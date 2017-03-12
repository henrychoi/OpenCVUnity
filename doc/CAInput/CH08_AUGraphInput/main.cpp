#include <stdio.h>
#include <AudioToolbox/AudioToolbox.h>
#include <Accelerate/Accelerate.h>
#include "CARingBuffer.h"

#pragma mark - utility functions -

// generic error handler - if err is nonzero, prints error message and exits program.
static void CheckError(OSStatus error, const char *operation)
{
    if (error == noErr) return;
    
    char str[20];
    // see if it appears to be a 4-char-code
    *(UInt32 *)(str + 1) = CFSwapInt32HostToBig(error);
    if (isprint(str[1]) && isprint(str[2]) && isprint(str[3]) && isprint(str[4])) {
        str[0] = str[5] = '\'';
        str[6] = '\0';
    } else
        // no, format it as an integer
        sprintf(str, "%d", (int)error);
    
    fprintf(stderr, "Error: %s (%s)\n", operation, str);
    
    exit(1);
}

typedef struct MyAUGraphPlayer
{
	AudioStreamBasicDescription streamFormat;
    AudioUnit inputUnit;
	AudioBufferList *inputBuffer;
    FILE* t_csv, *x_csv, *f_csv, *c_csv;
    float* padded_x;
    int8_t Eblock, iBlock;
} MyAUGraphPlayer;
MyAUGraphPlayer gPlayer = { .iBlock = -1 }
, *player = &gPlayer;

#define LOG2M 10
#define EXP_LOG2M (1<<LOG2M)
#define LOG2_HALFM (LOG2M - 1)
#define LOG2_2M (LOG2M + 1)

float FFT_pulse_conj_vDSPz_flattened[1<<LOG2_2M] __attribute__((aligned(16)))
= {
#include "FFT_pulse_conj_vDSPz_flattened.h"
};

static OSStatus InputRenderProc(void *inRefCon,
						 AudioUnitRenderActionFlags *ioActionFlags,
						 const AudioTimeStamp *inTimeStamp,
						 UInt32 inBusNumber,
						 UInt32 inNumberFrames,
						 AudioBufferList * ioData)
{
    if (player->iBlock < 0) { // The 1st block seems to have a pop sound for some reason
        player->iBlock++;
        return noErr;
    }
    
    //I only care about channel 0; copy out straight to the padded time
    //series so I can FFT
    player->inputBuffer->mBuffers[0].mData = (void*)
        (player->padded_x + (inNumberFrames * (2 * player->iBlock + 1)));

    CheckError(AudioUnitRender(player->inputUnit, ioActionFlags,
                            inTimeStamp, inBusNumber, inNumberFrames,
							player->inputBuffer),
               "AudioUnitRender failed");
    
    fprintf(player->t_csv, "%llu, %llu, %f, %u\n", // record the timestamps
            inTimeStamp->mHostTime, inTimeStamp->mWordClockTime,
            inTimeStamp->mRateScalar, inNumberFrames);
    if (++player->iBlock >= player->Eblock) {
        CheckError(AudioOutputUnitStop(player->inputUnit)
                   , "AudioOutputUnitStop");
    }
    return noErr;
}

#include <getopt.h>

static const char* optString = "Nctxfvh?";
static const struct option longOpts[] = {
    //{ "Filter", required_argument, NULL, 'F' },
    { "Nblock", required_argument, NULL, 'N' },
    { "Tcsv", required_argument, NULL, 't' },
    { "Xcsv", required_argument, NULL, 'x' },
    { "Fcsv", required_argument, NULL, 'f' },
    { "Ccsv", required_argument, NULL, 'c' },
    { "verbose", no_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};

static void die() {
    fprintf(stderr
            , "Usage: <-N Nblock> <-t time csv> <-x data csv> <-f f csv> <-c corr csv>\n"
            "\tNblock: number of callbacks to process\n"
            );
    exit(errno);
}

int main (int argc, char* const argv[]) {
    int longIndex = 0
        , opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'N':
                player->Eblock = (int8_t)atoi(optarg);
                break;
            case 'x':
                if (!(player->x_csv = fopen(optarg, "w"))) {
                    fprintf(stderr, "Failed to open %s", optarg);
                    die();
                }
                break;
            case 't':
                if (!(player->t_csv = fopen(optarg, "w"))) {
                    fprintf(stderr, "Failed to open %s", optarg);
                    die();
                }
                break;
            case 'f':
                if (!(player->f_csv = fopen(optarg, "w"))) {
                    fprintf(stderr, "Failed to open %s", optarg);
                    die();
                }
                break;
            case 'c':
                if (!(player->c_csv = fopen(optarg, "w"))) {
                    fprintf(stderr, "Failed to open %s", optarg);
                    die();
                }
                break;
            default:
                die();
        }
        
        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }
    
    // generate description that will match audio HAL
    AudioComponentDescription inputcd = {0};
    inputcd.componentType = kAudioUnitType_Output;
    inputcd.componentSubType = kAudioUnitSubType_HALOutput;
    inputcd.componentManufacturer = kAudioUnitManufacturer_Apple;
    
    AudioComponent comp = AudioComponentFindNext(NULL, &inputcd);
    if (comp == NULL) {
        printf ("can't get output unit");
        exit (-1);
    }
    
    CheckError(AudioComponentInstanceNew(comp, &player->inputUnit),
               "Couldn't open component for inputUnit");
    
    // enable/io
    UInt32 disableFlag = 0;
    UInt32 enableFlag = 1;
    AudioUnitScope outputBus = 0;
    AudioUnitScope inputBus = 1;
    CheckError (AudioUnitSetProperty(player->inputUnit,
                                     kAudioOutputUnitProperty_EnableIO,
                                     kAudioUnitScope_Input,
                                     inputBus,
                                     &enableFlag,
                                     sizeof(enableFlag)),
                "Couldn't enable input on I/O unit");
    
    CheckError (AudioUnitSetProperty(player->inputUnit,
                                     kAudioOutputUnitProperty_EnableIO,
                                     kAudioUnitScope_Output,
                                     outputBus,
                                     &disableFlag,	// well crap, have to disable
                                     sizeof(enableFlag)),
                "Couldn't disable output on I/O unit");
    
    // set device (osx only... iphone has only one device)
    AudioDeviceID defaultDevice = kAudioObjectUnknown;
    UInt32 propertySize = sizeof (defaultDevice);
    
    // AudioHardwareGetProperty() is deprecated
    //	CheckError (AudioHardwareGetProperty(kAudioHardwarePropertyDefaultInputDevice,
    //										 &propertySize,
    //										 &defaultDevice),
    //				"Couldn't get default input device");
    
    // AudioObjectProperty stuff new in 10.6, replaces AudioHardwareGetProperty() call
    // TODO: need to update ch08 to explain, use this call. need CoreAudio.framework
    AudioObjectPropertyAddress defaultDeviceProperty;
    defaultDeviceProperty.mSelector = kAudioHardwarePropertyDefaultInputDevice;
    defaultDeviceProperty.mScope = kAudioObjectPropertyScopeGlobal;
    defaultDeviceProperty.mElement = kAudioObjectPropertyElementMaster;
    
    CheckError (AudioObjectGetPropertyData(kAudioObjectSystemObject,
                                           &defaultDeviceProperty,
                                           0,
                                           NULL,
                                           &propertySize,
                                           &defaultDevice),
                "Couldn't get default input device");
    
    // set this defaultDevice as the current device
    // kAudioUnitErr_InvalidPropertyValue if output is enabled on inputUnit
    CheckError(AudioUnitSetProperty(player->inputUnit,
                                    kAudioOutputUnitProperty_CurrentDevice,
                                    kAudioUnitScope_Global,
                                    outputBus,
                                    &defaultDevice,
                                    sizeof(defaultDevice)),
               "Couldn't set default device as the current (output unit) device");
    
    propertySize = sizeof (AudioStreamBasicDescription);
    CheckError(AudioUnitGetProperty(player->inputUnit,
                                    kAudioUnitProperty_StreamFormat,
                                    kAudioUnitScope_Output, //to other units
                                    inputBus,
                                    &player->streamFormat,
                                    &propertySize),
               "Couldn't get ASBD from output unit");
    //player->streamFormat.mChannelsPerFrame = 1; // want mono
    //Filter designed for this rate, so fix at expected frequency
    player->streamFormat.mSampleRate = 44100;
    CheckError(AudioUnitSetProperty(player->inputUnit,
                                    kAudioUnitProperty_StreamFormat,
                                    kAudioUnitScope_Output,
                                    inputBus,
                                    &player->streamFormat,
                                    propertySize),
               "Couldn't set output unit Fs");
    
    // 9/6/10 - check the input device's stream format
    AudioStreamBasicDescription deviceFormat;
    CheckError(AudioUnitGetProperty(player->inputUnit,
                                    kAudioUnitProperty_StreamFormat,
                                    kAudioUnitScope_Input, //from HW to IO unit
                                    inputBus,
                                    &deviceFormat,
                                    &propertySize),
               "Couldn't get ASBD from input unit");
    
    deviceFormat.mSampleRate = 44100;
    CheckError(AudioUnitSetProperty(player->inputUnit,
                                    kAudioUnitProperty_StreamFormat,
                                    kAudioUnitScope_Input,
                                    inputBus,
                                    &deviceFormat,
                                    propertySize),
               "Couldn't set Fs input unit");
    
    /* allocate some buffers to hold samples between input and output callbacks
     (this part largely copied from CAPlayThrough) */
    //Get the size of the IO buffer(s)
    UInt32 bufferSizeFrames = 1024;
    propertySize = sizeof(UInt32);
    CheckError(AudioUnitSetProperty(player->inputUnit,
                                    kAudioDevicePropertyBufferFrameSize,
                                    kAudioUnitScope_Global,
                                    0,
                                    &bufferSizeFrames,
                                    propertySize),
               "Couldn't set buffer frame size to input unit");

    CheckError (AudioUnitGetProperty(player->inputUnit,
                                     kAudioDevicePropertyBufferFrameSize,
                                     kAudioUnitScope_Global,
                                     0,
                                     &bufferSizeFrames,
                                     &propertySize),
                "Couldn't get buffer frame size from input unit");
    if (bufferSizeFrames != 1024) {
        fprintf(stderr, "bufferSizeFrames change didn't stick");
        exit(1);
    }
    
    UInt32 bufferSizeBytes = bufferSizeFrames * sizeof(Float32);
    //The time series data for channel 0 should wind up here.
    player->padded_x = (float*)calloc(2 * player->Eblock, bufferSizeBytes);
    // I think the CoreAudio book is wrong about the buffer size in case of
    // interleaved format.  Just supporting the non-interleaved case until I
    // understand the interleaved format
    if (!(player->streamFormat.mFormatFlags & kAudioFormatFlagIsNonInterleaved)) {
        fprintf(stderr, "not what I wanted: format is interleaved\n");
        return false;
    }
    // allocate an AudioBufferList plus enough space for array of AudioBuffers
    UInt32 propsize = offsetof(AudioBufferList, mBuffers[0])
    + (sizeof(AudioBuffer) * player->streamFormat.mChannelsPerFrame);
    
    //malloc buffer lists
    player->inputBuffer = (AudioBufferList *)malloc(propsize);
    player->inputBuffer->mNumberBuffers = player->streamFormat.mChannelsPerFrame;
    //pre-malloc buffers for AudioBufferLists
    UInt32 i=0;
    player->inputBuffer->mBuffers[i].mNumberChannels = 1;
    player->inputBuffer->mBuffers[i].mDataByteSize = bufferSizeBytes;
    for(i=1; i< player->inputBuffer->mNumberBuffers; i++) {
        player->inputBuffer->mBuffers[i].mNumberChannels = 1;
        player->inputBuffer->mBuffers[i].mDataByteSize = bufferSizeBytes;
        player->inputBuffer->mBuffers[i].mData = malloc(bufferSizeBytes);
    }
    
    // set render proc to supply samples from input unit
    AURenderCallbackStruct callbackStruct;
    callbackStruct.inputProc = InputRenderProc; 
    callbackStruct.inputProcRefCon = player;
    
    CheckError(AudioUnitSetProperty(player->inputUnit, 
                                    kAudioOutputUnitProperty_SetInputCallback, 
                                    kAudioUnitScope_Global,
                                    0,
                                    &callbackStruct, 
                                    sizeof(callbackStruct)),
               "Couldn't set input callback");
    
    CheckError(AudioUnitInitialize(player->inputUnit),
               "Couldn't initialize input unit");
    // Input unit is now ready for use
    
    CheckError(AudioOutputUnitStart(player->inputUnit)
               , "AudioOutputUnitStart");

    //Normally, time series of length M only requires a split buffer of M/2 for the
    //real and complex portions (so that the total number is still M).  But because I
    //zero pad, I need double this amount.  For M=1024, I will eat up 16 KB for these
    //temporary vars, which would not do on an MCU.
    float creal[1<<LOG2M] __attribute__((aligned(16))) //input/output buffer
    , cimag[1<<LOG2M] __attribute__((aligned(16)))
    , treal[1<<LOG2M] __attribute__((aligned(16))) //temporary scratch pad
    , timag[1<<LOG2M] __attribute__((aligned(16)))
    ;
    COMPLEX_SPLIT splitc = { .realp = creal, .imagp = cimag }
    , splitt = { .realp = treal, .imagp = timag }
    , split_filter = { .realp = FFT_pulse_conj_vDSPz_flattened
        , .imagp = FFT_pulse_conj_vDSPz_flattened + (1<<LOG2M) }
    ;
#if 0
    // Write the filter coefficients for sanity check
    for(UInt32 k=0; k < (1<<LOG2M); ++k) { // Output to CSV for Matlab viewing
        fprintf(player->f_csv, "%f,%f\n", split_filter.realp[k], split_filter.imagp[k]);
    }
#endif
    //save the filter Nyquist freq response
    float FFT_filter_nyq = *split_filter.imagp; *split_filter.imagp = 0;

    for (int8_t iBlock = 0; iBlock < player->Eblock; ) {
        if (iBlock >= player->iBlock) {
            usleep(1000);
            continue;
        }
        // Have a new block to process
        float* padded_x = player->padded_x + 2 * bufferSizeFrames * iBlock
            , *x = padded_x;

        for(UInt32 k=0; k < 2 * bufferSizeFrames; ++k)
            fprintf(player->x_csv, "%f\n", *x++);

        vDSP_ctoz((COMPLEX*)padded_x, 2, &splitc, 1, 1<<LOG2M);
        
        static FFTSetup fftSetup = vDSP_create_fftsetup(LOG2_2M, kFFTRadix2);
        vDSP_fft_zript(fftSetup, &splitc, 1, &splitt, LOG2_2M, FFT_FORWARD);
        // splitc now has FFT(x)--but in vDSP packed real format (only +f spectrum)
        // The FFT(x)[M], which is purely real, is in the imag(FFT(x)[0].  The same
        // goes for split_filter.  I must save those real values and set the
        // img(FFT(x)[k=0]) to 0 for correct result.
        float FFT_c_nyq = *splitc.imagp; *splitc.imagp = 0;
#if 1
        vDSP_zvmul(&splitc, 1, &split_filter, 1, &splitc, 1
                   , 1<<LOG2M//The multiplication is over vector of length M
                   , 1 // don't conjugate (even though we are correlating);
                   ); // because the filter is already in conjugate form
        // Restore the Nyquist frequency portion, which is pure real
        *splitc.imagp = FFT_c_nyq * FFT_filter_nyq; //Restore the Nyquist frequency portion
#endif
        for(UInt32 k=0; k < (1<<LOG2M); ++k) { // Output to CSV for Matlab viewing
            fprintf(player->f_csv, "%f,%f\n", splitc.realp[k], splitc.imagp[k]);
        }
        vDSP_fft_zript(fftSetup, &splitc, 1, &splitt, LOG2_2M, FFT_INVERSE);
        //vDSP_fft_zrip(fftSetup, &splitc, 1, LOG2_2M, FFT_INVERSE);
        vDSP_ztoc(&splitc, 1, (COMPLEX*)padded_x, 2, 1<<LOG2M);// Go back to real space
        //padded_x[0] = 1E9; //debug the line marker
        for(UInt32 k=0; k < (1<<LOG2_2M); ++k) { // Output to CSV for Matlab viewing
            fprintf(player->c_csv, "%f\n", padded_x[k]);
        }
        
        iBlock++;
    }
cleanup:
    fclose(player->c_csv);
    fclose(player->t_csv);
    fclose(player->x_csv);
    fclose(player->f_csv);
    return 0;
}
