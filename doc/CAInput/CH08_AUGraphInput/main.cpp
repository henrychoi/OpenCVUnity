#include <stdio.h>
#include <AudioToolbox/AudioToolbox.h>
#include <Accelerate/Accelerate.h>
#include "llsMP.h"
#include "llsQ.h"

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

typedef struct MyMic
{
	AudioStreamBasicDescription streamFormat;
    AudioUnit inputUnit;
	AudioBufferList *inputBuffer;
    FILE* t_csv, *x_csv, *f_csv, *c_csv;
    struct llsMP mpool;
    struct llsQ padded_xQ;
    int8_t Eblock, iBlock;
    bool firstBlock;
} MyMic;
MyMic gPlayer = { .firstBlock = true }
, *player = &gPlayer;

#define LOG2M 10
#define EXP_LOG2M (1<<LOG2M)
#define LOG2_HALFM (LOG2M - 1)
#define LOG2_2M (LOG2M + 1)

// The matched filter from Matlab
float FFT_pulse_conj_vDSPz_flattened[1<<LOG2_2M] __attribute__((aligned(16)))
= {
#include "FFT_pulse_conj_vDSPz_flattened.h"
};

struct SampleBlock {
    UInt64 mHostTime;//The timestamp for these frames
    //UInt64 mWordClockTime //don't bother; always zero

    //I can calculate the time of the 1st sample with this just as well
    uint32_t iBlock;
    UInt32 nFrames;
    //float mRateScalar; //had 4 bytes left, so record the coarse version of scalar
    float sample[1<<LOG2_2M];
};

static OSStatus InputRenderProc(void *inRefCon,
						 AudioUnitRenderActionFlags *ioActionFlags,
						 const AudioTimeStamp *inTimeStamp,
						 UInt32 inBusNumber,
						 UInt32 inNumberFrames,
						 AudioBufferList * ioData)
{
    if (player->firstBlock) { // The 1st block seems to have a pop sound,
        // at least on my iPhone 6s.  Throw it away.  There will be plenty
        // more blocks.
        player->firstBlock = false;
        return noErr;
    }
    
    //I only care about channel 0 (Q: is that L or R?  I think L); copy it
    //out straight to the padded time series so I can FFT
    struct SampleBlock* block = (struct SampleBlock*)llsMP_get(&player->mpool);
    if (!block) {
        fprintf(stderr, "Memory pool exhausted\n");
        return errno;
    }
    player->inputBuffer->mBuffers[0].mData = (void*)(&block->sample[1<<LOG2M]);

    CheckError(AudioUnitRender(player->inputUnit, ioActionFlags,
                            inTimeStamp, inBusNumber, inNumberFrames,
							player->inputBuffer),
               "AudioUnitRender failed");

    block->iBlock = player->iBlock;
    block->mHostTime = inTimeStamp->mHostTime;
    block->nFrames = inNumberFrames;
    if (!llsQ_push(&player->padded_xQ, block)) {
        fprintf(stderr, "llsQ_push failed\n");
        return errno;
    }
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
    uint32_t k;

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
        fprintf (stderr, "can't get output unit");
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
    UInt32 bufferSizeFrames = 1<<LOG2M;
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
    if (bufferSizeFrames != 1<<LOG2M) {
        fprintf(stderr, "bufferSizeFrames change didn't stick");
        exit(1);
    }
    
    UInt32 bufferSizeBytes = bufferSizeFrames * sizeof(Float32);
    // Memory pool to hold the padded x
    if (!llsMP_alloc(&player->mpool
                     , 3 // capacity
                     , sizeof(struct SampleBlock) // memsize
                     , 16 // alignment: recommended by vDSP programming guide
                     , 1)//memset zero, since this pool is used for zero padding
        ) {
        fprintf (stderr, "Can't allocate zero padded memory pool");
        exit (errno);
    }
    
#undef UNIT_TEST_LLSMP
#ifdef UNIT_TEST_LLSMP
    for (k=0; k < 3; ++k) {
        struct SampleBlock* block = (struct SampleBlock*)llsMP_get(&player->mpool);
        if (!block) {
            fprintf(stderr, "Memory pool exhausted\n");
            return errno;
        }
        float* x = block->sample;
        for(uint32_t j=0; j < (1<<LOG2_2M); ++j)
            fprintf(player->x_csv, "%f\n", *x++);
    }
    return noErr;
#endif
    
    if (!llsQ_alloc(&player->padded_xQ, 1)) {
        fprintf (stderr, "Can't allocate memory Q");
        exit (errno);
    }

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
        // I won't even look at the other channels, but the Render method
        // still wants some place to copy the samples to
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
        , ftemp[1<<LOG2_2M] __attribute__((aligned(16))) //temporary scratch pad
        , overlap_save[1<<LOG2M] __attribute__((aligned(16))) //final correlation
    ;
    memset(overlap_save, 0, sizeof(overlap_save));
    COMPLEX_SPLIT splitc = { .realp = creal, .imagp = cimag }
        , split_temp = { .realp = ftemp, .imagp = &ftemp[1<<LOG2M] }
        , split_filter = { .realp = FFT_pulse_conj_vDSPz_flattened
            , .imagp = FFT_pulse_conj_vDSPz_flattened + (1<<LOG2M) }
    ;
#if 0
    // Write the filter coefficients for sanity check
    for(k=0; k < (1<<LOG2M); ++k) { // Output to CSV for Matlab viewing
        fprintf(player->f_csv, "%f,%f\n", split_filter.realp[k], split_filter.imagp[k]);
    }
#endif
    //save the filter Nyquist freq response
    float FFT_filter_nyq = *split_filter.imagp; *split_filter.imagp = 0;

    for (int8_t iBlock = 0; iBlock < player->Eblock; ) {
        struct SampleBlock* block;
        if (!llsQ_pop(&player->padded_xQ, (void**)&block)) {
            usleep(1000);
            continue;
        }
        // Have a new block to process
        fprintf(player->t_csv, "%u, %llu\n", // record the timestamps
                block->iBlock, block->mHostTime);

        for(float* x = &block->sample[1<<LOG2M], k=0; k < (1<<LOG2M); ++k)
            fprintf(player->x_csv, "%f\n", *x++);

        vDSP_ctoz((COMPLEX*)block->sample, 2, &splitc, 1, 1<<LOG2M);

        // Copied out the sample, so I can return
        if (!llsMP_return(&player->mpool, block)) {// return the memory to mpool
            fprintf(stderr, "llsMP_return failed");
            exit(errno);
        }
        
        // Rest of the calculation can be done with stack
        static FFTSetup fftSetup = vDSP_create_fftsetup(LOG2_2M, kFFTRadix2);
        vDSP_fft_zript(fftSetup, &splitc, 1, &split_temp, LOG2_2M, FFT_FORWARD);
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
        for(k=0; k < (1<<LOG2M); ++k) { // Output to CSV for Matlab viewing
            fprintf(player->f_csv, "%f,%f\n", splitc.realp[k], splitc.imagp[k]);
        }
        vDSP_fft_zript(fftSetup, &splitc, 1, &split_temp, LOG2_2M, FFT_INVERSE);
        //vDSP_fft_zrip(fftSetup, &splitc, 1, LOG2_2M, FFT_INVERSE);
        
        // Go back to real space, but I must not write back into the zero
        // padded sample buffer: this will ruin zero pading!
        vDSP_ztoc(&splitc, 1, (COMPLEX*)ftemp, 2, 1<<LOG2M);
        
        //Correlation to output: overlap_save + ftemp[0:M-1];
        vDSP_vadd(ftemp, 1, overlap_save, 1, ftemp, 1, 1<<LOG2M);
        
        //Correlation to save for next round: overlap_save = ftemp[M:2M-1]
        memcpy(overlap_save, &ftemp[1<<LOG2M], sizeof(overlap_save));
        
        //padded_x[0] = 1E9; //debug the line marker
        for(k=0; k < (1<<LOG2M); ++k) { // Output correlation to CSV for Matlab viewing
            fprintf(player->c_csv, "%f\n", ftemp[k]);
        }
        ++iBlock;
    }
cleanup:
    llsMP_free(&player->mpool);
    fclose(player->c_csv);
    fclose(player->t_csv);
    fclose(player->x_csv);
    fclose(player->f_csv);
    return 0;
}
