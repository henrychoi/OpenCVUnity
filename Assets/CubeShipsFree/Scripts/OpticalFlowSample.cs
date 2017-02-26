using UnityEngine;
using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
//using System.Diagnostics;
using CubeSpaceFree;//coupled to the 3D_Example/PlayerController
using UnityEngine.Assertions;

#if UNITY_5_3 || UNITY_5_3_OR_NEWER
using UnityEngine.SceneManagement;
#endif
using OpenCVForUnity;

namespace OpenCVForUnitySample
{
    /// <summary>
    /// Optical flow sample.
    /// http://stackoverflow.com/questions/6505779/android-optical-flow-with-opencv?rq=1
    /// </summary>
    public class OpticalFlowSample : MonoBehaviour
    {
        
        /// <summary>
        /// The mat op flow this.
        /// </summary>
        Mat matOpFlowThis;

        /// <summary>
        /// The mat op flow previous.
        /// </summary>
        Mat matOpFlowPrev;

        /// <summary>
        /// The i GFFT max.
        /// </summary>
        int iGFFTMax = 40;

        /// <summary>
        /// The MO pcorners.
        /// </summary>
        MatOfPoint MOPcorners;

        /// <summary>
        /// The m MO p2fpts this.
        /// </summary>
        MatOfPoint2f mMOP2fptsThis;

        /// <summary>
        /// The m MO p2fpts previous.
        /// </summary>
        MatOfPoint2f mMOP2fptsPrev;

        /// <summary>
        /// The m MO p2fpts safe.
        /// </summary>
        MatOfPoint2f mMOP2fptsSafe;

        /// <summary>
        /// The m MOB status.
        /// </summary>
        MatOfByte mMOBStatus;

        /// <summary>
        /// The m MO ferr.
        /// </summary>
        MatOfFloat mMOFerr;

        /// <summary>
        /// The color red.
        /// </summary>
        Scalar colorRed = new Scalar (255, 0, 0, 255);

        /// <summary>
        /// The i line thickness.
        /// </summary>
        int iLineThickness = 3;

        /// <summary>
        /// The texture.
        /// </summary>
        Texture2D texture;

        /// <summary>
        /// The web cam texture to mat helper.
        /// </summary>
        WebCamTextureToMatHelper webCamTextureToMatHelper;

		PlayerController playerCtrl;

		const string camCalibrationFn = "camera_intrinsic.txt";
		string camCalFilepath;
    
        // Use this for initialization
        void Start ()
        {
            webCamTextureToMatHelper = gameObject.GetComponent<WebCamTextureToMatHelper> ();
            webCamTextureToMatHelper.Init ();
    
			//gameObject is what this script is attached to.  In this project, it is the Quad
			//Quad's parent is the Player, which has the PlayerController script component
			playerCtrl = gameObject.GetComponentInParent<PlayerController> ();
			Assert.IsNotNull (playerCtrl);

			//Debug.Log ("Found yml " + blobparams_yml_filepath);
			blobDetector = FeatureDetector.create (FeatureDetector.SIMPLEBLOB);
			blobDetector.read (Utils.getFilePath ("blobparams.yml"));

			//Read stored camera calibration /////////////////////////////////
			camCalFilepath = Utils.getFilePath (camCalibrationFn);
			bool cameraCalParseSuccess = true;

			if (!File.Exists(camCalFilepath)) {
				cameraCalParseSuccess = false;
				goto doneCameraCalParse;
			}
				
			string[] lines = File.ReadAllLines(camCalFilepath);
			Debug.Log ("Read " + lines.Length + " lines from " + camCalFilepath);
			if (lines.Length < (2+3*3 + 2+5)) {//intrinsic and distortion should be separate lines
				Debug.Log("Found only " + lines.Length + " lines in " + camCalFilepath);
				cameraCalParseSuccess = false;
				goto doneCameraCalParse;
			}
			int typ;
			int iLine=0;
			//Skip the 1st line (only a comment)
			if (!Int32.TryParse(lines[++iLine], out typ)) {
				Debug.Log("Failed to parse line " + iLine + ": " + lines[iLine]);
				cameraCalParseSuccess = false;
				goto doneCameraCalParse;
			}
			intrinsic_matrix = new Mat(3, 3, typ);
			for(int i=0; i < 3; ++i) {
				for(int j=0; j < 3; ++j) {
					double d;
					if (!Double.TryParse(lines[++iLine], out d)) {
						Debug.Log("Failed to parse line " + iLine + ": " + lines[iLine]);
						cameraCalParseSuccess = false;
						goto doneCameraCalParse;
					}
					intrinsic_matrix.put(i, j, d);
				}
			}
			++iLine;// Skip 1 line before distortion mx
			if (!Int32.TryParse(lines[++iLine], out typ)) {
				Debug.Log("Failed to parse line " + iLine + ": " + lines[iLine]);
				cameraCalParseSuccess = false;
				goto doneCameraCalParse;
			}
			distortion_coeffs = new Mat(1, 5, typ);
			for(int j=0; j < 5; ++j) {
				double d;
				if (!Double.TryParse(lines[++iLine], out d)) {
					Debug.Log("Failed to parse line " + iLine + ": " + lines[iLine]);
					cameraCalParseSuccess = false;
					goto doneCameraCalParse;
				}
				distortion_coeffs.put(0, j, d);
			}
		doneCameraCalParse:
			if (cameraCalParseSuccess) { // show the read mx
				Debug.Log (String.Format ("{0} -> intrinsic {1}, distortion {2}",
					camCalFilepath,
					intrinsic_matrix.dump(), //came back as 3x3
					distortion_coeffs.dump()) //came back as 1x5
				);
			} else {
				intrinsic_matrix = new Mat(); //invalidate the intrinsic mx
			} //////////////////////////////////////////////////////////////////////////

			// Dead reckon the physical corners of the chessboard
			const float chessboard_scale = 0.026f; // Each chessboard square is a 26 mm long
			List<Point3> corners = new List<Point3>(7*7);
			 //Dead reckon the chessboard coner locations
			for (int i = 0; i < 7; ++i) {
				for (int j = 0; j < 7; ++j) {
					corners.Add(new Point3(i*chessboard_scale, j*chessboard_scale, 0));
				}
			}
			MatOfPoint3f wld_corners = new MatOfPoint3f();
			wld_corners.fromList (corners);
			for (int i = 0; i < nBoard; ++i) { // Will use the same board for all 10 images
				wld_corner_list.Add (wld_corners);
			}
			Debug.Log (String.Format ("Added {0} wld_corners", wld_corner_list.Count));
        }

		const float iPhone6s_FoV = Mathf.Deg2Rad * 65/2;

        /// <summary>
        /// Raises the web cam texture to mat helper inited event.
        /// </summary>
        public void OnWebCamTextureToMatHelperInited ()
        {
            Debug.Log ("OnWebCamTextureToMatHelperInited");
            
            Mat webCamTextureMat = webCamTextureToMatHelper.GetMat ();

            texture = new Texture2D (webCamTextureMat.cols (), webCamTextureMat.rows (), TextureFormat.RGBA32, false);

			//gameObject is what this script is attached to.  In this project, it is the Quad
			//gameObject.GetComponent<Renderer> is the Renderer that is attached to the gameObject:
			//	MeshRenderer in this project
			//gameObject.GetComponent<Renderer>.material is the *1st* material attached to the renderer;
			//	In this project, the MeshRenderer only has 1 material in its materials list: the Missing(Material) 
            gameObject.GetComponent<Renderer> ().material.mainTexture = texture;

			//localScale is the scale of the transform relative to the parent.  The webcamTexture is sized by the
			//image returned from the camera, for example 1280x720, 1920x1080.  I need to fill the whole surface
			//of the Quad sitting 2000 [screen units] away from the main perspective camera with a FoV = 65 deg
			//(and a far clipping plane at 1000).  The FoV may be programmatically queriable, in which case the
			//equation will change.  From a simple right triangle with the bottom leg = 2000 and the angle FoV/2,
			//the desired width of the image is 2000 * 2 * tan(FoV/2), which is ~2500.  The scale factor to map
			//the webCamTexture height to te Quad height is 2000 * tan(FoV/2) / webCamTextureMat.rows ()
			float d2quad = gameObject.transform.position.magnitude;//player orientation is all messed up, so use magnitude
			//Debug.Log("Quad transform @" + d2quad);
			Assert.IsTrue(d2quad > 1000);
			float texture2quadHeight = d2quad * 2 * Mathf.Tan(iPhone6s_FoV);
			int height = webCamTextureMat.height (), width = webCamTextureMat.width ();
			gameObject.transform.localScale = new Vector3((texture2quadHeight / height) * width, texture2quadHeight, 1);
			//Debug.Log ("localScale = " + gameObject.transform.localScale);
            
            /* Following block is unnecessary if perspective projection is used for the camera
			Debug.Log ("Screen.width " + Screen.width + " Screen.height " + Screen.height
				+ " Screen.orientation " + Screen.orientation);

            float width = webCamTextureMat.width();
            float height = webCamTextureMat.height();
            
            float widthScale = (float)Screen.width / width;
            float heightScale = (float)Screen.height / height;
            if (widthScale < heightScale) {
                Camera.main.orthographicSize = (width * (float)Screen.height / (float)Screen.width) / 2;
            } else {
                Camera.main.orthographicSize = height / 2;
            }
			*/
			//Mat (height, width, int type, Scalar s)
            matOpFlowThis = new Mat ();
            matOpFlowPrev = new Mat ();
            MOPcorners = new MatOfPoint ();
            mMOP2fptsThis = new MatOfPoint2f ();
            mMOP2fptsPrev = new MatOfPoint2f ();
            mMOP2fptsSafe = new MatOfPoint2f ();
            mMOBStatus = new MatOfByte ();
            mMOFerr = new MatOfFloat ();
        }

        /// <summary>
        /// Raises the web cam texture to mat helper disposed event.
        /// </summary>
        public void OnWebCamTextureToMatHelperDisposed ()
        {
            Debug.Log ("OnWebCamTextureToMatHelperDisposed");

			//if (rxgxbMat != null) rxgxbMat.Dispose ();
			if (tempA != null) tempA.Dispose ();
			if (tempB != null) tempB.Dispose ();

			if (matOpFlowThis != null) matOpFlowThis.Dispose ();
            if (matOpFlowPrev != null) matOpFlowPrev.Dispose ();
            if (MOPcorners != null) MOPcorners.Dispose ();
            if (mMOP2fptsThis != null) mMOP2fptsThis.Dispose ();
            if (mMOP2fptsPrev != null) mMOP2fptsPrev.Dispose ();
            if (mMOP2fptsSafe != null) mMOP2fptsSafe.Dispose ();
            if (mMOBStatus != null) mMOBStatus.Dispose ();
            if (mMOFerr != null) mMOFerr.Dispose ();
        }

        /// <summary>
        /// Raises the web cam texture to mat helper error occurred event.
        /// </summary>
        /// <param name="errorCode">Error code.</param>
        public void OnWebCamTextureToMatHelperErrorOccurred(WebCamTextureToMatHelper.ErrorCode errorCode){
            Debug.Log ("OnWebCamTextureToMatHelperErrorOccurred " + errorCode);
        }

		float lastOdoTime = .0f, lastChessboardTime = .0f;

		Mat tempA, tempB,// rxgxbMat,
			emptyMat = new Mat()
			//erode_kernel = new Mat()
				//Imgproc.getStructuringElement (Imgproc.MORPH_ELLIPSE, new Size(7,7))
		;
		Point erode_anchor = new Point(0,0);

		FeatureDetector blobDetector;
		Scalar channelMultiplyScale = new Scalar(1.0f/256);
		Size chessboard_size = new Size(7, 7);//standard OpenCV square chessboard
		//Q: why can't I make this a const?

		const int nBoard = 10; //Look for 10 successful chessboard images with all corners found
		List<Mat> img_corner_list = new List<Mat>(nBoard);
		List<Mat> wld_corner_list = new List<Mat>(nBoard);//pre-populated in Start()
		Mat intrinsic_matrix = new Mat(),
		    distortion_coeffs = new Mat();
		System.Diagnostics.Stopwatch watch = new System.Diagnostics.Stopwatch();
        // Update is called once per frame
        void Update ()
        {
			if (!(webCamTextureToMatHelper.IsPlaying () && webCamTextureToMatHelper.DidUpdateThisFrame ()))
				return;

			float now = Time.time, dT = now - lastOdoTime;
			if (dT < 0.1f) // Don't run visual odometry/sensor-fusion too fast (sucks current)
				return;
			lastOdoTime = now;

			Mat rgbaMat = webCamTextureToMatHelper.GetMat ();
			int width = rgbaMat.width (), height = rgbaMat.height ();
			//watch.Start ();
			// Assume the brightest white region(s) are the torches
			List<Mat> channels = new List<Mat>();
			Core.split (rgbaMat, channels);
			//Debug.Log ("split channel type = " + channels [0].type ());

			if (tempA == null) tempA = new Mat (height, width, channels [0].type ());
			if (tempB == null) tempB = new Mat (height, width, channels [0].type ());
			//if (rxgxbMat == null) rxgxbMat = new Mat (height, width, channels [0].type ());

			//Debug.Log ("channels: " + channels.Count);
			//Core.mulSpectrums(channels[0], channels[1], rxgxbMat
			Mat rxgxbMat = channels[0].mul(channels[1].mul(channels[2], 1.0f/256), 1.0f/256);
			//Core.multiply (channels[1], channels[2], tempA, 1.0f/256);
			//Core.multiply (tempA, channels [0], tempB, 1.0f/256);

			Imgproc.medianBlur (rxgxbMat, tempA, 7);
			Imgproc.morphologyEx (tempA, tempB, Imgproc.MORPH_ERODE, emptyMat, erode_anchor, 2);
			Core.MinMaxLocResult minmax = Core.minMaxLoc (tempB);
			//Debug.Log (String.Format("max {0} @ {1}", minmax.maxVal, minmax.maxLoc));

			//Wipe pixels significantly dimmer than the max.  Initiailly I wanted to
			//change the blob detector's minThreshold dynamically, but since the params are in
			//a file, it's too onerous to change for every Update.
			Core.subtract (tempB,
				new Scalar(Mathf.Max((float)(0.9f * minmax.maxVal), 200.0f)),
				tempA);
			//Mat normalized = tempB.mul((256.0f/0.25f) * minmax.maxVal);
			//Core.multiply(tempA, new Scalar((128.0f/0.2f) * minmax.maxVal), tempB);

			MatOfKeyPoint keypoints = new MatOfKeyPoint ();
			blobDetector.detect (tempA, keypoints);
			//Debug.Log ("keypoints found " + keypoints.size ());
			//watch.Stop ();
			//Debug.Log ("blob detection took " + watch.Elapsed);
			//watch.Reset ();
			/*
			Imgproc.calcHist (new List<Mat> (new Mat[]{ rxgxbMat }),
			    new MatOfInt (0), emptyMat,empty maskMat => no mask
				histMat, //output
				new MatOfInt (16), //histSize(s): # bins
				new MatOfFloat (0, 180)//the left and right most bins
			);
			Debug.Log ("hist " + histMat);
			*/

			if (now - lastChessboardTime < 1.0f)
				goto afterChessboard;
			
			// Find chessboard
			watch.Start ();
			MatOfPoint2f corners = new MatOfPoint2f();
			bool found = Calib3d.findChessboardCorners (rgbaMat, chessboard_size, corners,
				Calib3d.CALIB_CB_ADAPTIVE_THRESH | Calib3d.CALIB_CB_NORMALIZE_IMAGE //OpenCV default
				| Calib3d.CALIB_CB_FAST_CHECK);//no guarantee that I will hold up the chessboard, so skip fast
			watch.Stop ();
			Calib3d.drawChessboardCorners(rgbaMat, chessboard_size, corners, found); //Let's see those corners!
			Debug.Log ("chessboard detection took " + watch.Elapsed);
			watch.Reset ();

			if (found) {
				//Debug.Log (String.Format("Found {0}/{1} corners", img_corner_list.Count, nBoard));
				if (intrinsic_matrix.dims() != 0) //Already have the camera calibration
					goto chessboardPnP;

				// Else need to calibrate
				img_corner_list.Add(corners);
				if (img_corner_list.Count < nBoard)
					goto afterChessboard;
				// Have enough imag coners for calibration
				List<Mat> rvecs = new List<Mat>(), tvecs = new List<Mat>();
				double camCalError = Calib3d.calibrateCamera (
					wld_corner_list, img_corner_list, rgbaMat.size (), //input
					intrinsic_matrix, distortion_coeffs, rvecs, tvecs); //output
				Debug.Log (String.Format (
					"calibrate error {0}\nintrinsic type {1} mx {2}\ndistortion type {3} mx {4}",
					camCalError, //on the order of 1 pixel
					intrinsic_matrix.type(), intrinsic_matrix.dump(), //came back as 3x3
					distortion_coeffs.type(), distortion_coeffs.dump()) //came back as 1x5
				);

				/* save the intrinsic params--to where?  StreamingAssets doesn't seem to work
				Debug.Log("Writing calibration file " + camCalFilepath);
				StreamWriter writer = File.CreateText(camCalFilepath);

				writer.Write (intrinsic_matrix.type ()); writer.WriteLine ();
				for (int i = 0; i < 3; ++i) {
					double[] row = intrinsic_matrix.get (i, 0);
					Debug.Log ("Intrinsic row " + row);
					for (int j = 0; j < 3; ++j) {
						writer.Write (row [j]); writer.WriteLine ();
					}
				}
				writer.WriteLine (); // a separator line

				writer.Write (distortion_coeffs.type ()); writer.WriteLine ();
				double[] d = distortion_coeffs.get (0, 0);
				Debug.Log ("distortion " + d);
				for (int j = 0; j < 5; ++j) {
					writer.Write (d [j]); writer.WriteLine ();
				}
				writer.Close ();
				*/

			chessboardPnP:
				/*
				 Calib3d.solvePnP	(	MatOfPoint3f 	objectPoints, MatOfPoint2f 	imagePoints,
					Mat 	cameraMatrix,
					MatOfDouble 	distCoeffs,
					Mat 	rvec,
					Mat 	tvec 
					)	
				 */
				;
			}
			lastChessboardTime = now;
	afterChessboard:
			//Debug.Log("max of cross image = " + rxgxbMat.
			//Debug.Log (String.Format ("original size {0}x{1}", rgbaMat.width(), rgbaMat.height()));
			//Debug.Log (String.Format ("new size {0}x{1}", rxgxbMat.width(), rxgxbMat.height()));

			Features2d.drawKeypoints (rgbaMat, keypoints, rgbaMat);

			Utils.matToTexture2D (rgbaMat, texture, webCamTextureToMatHelper.GetBufferColors());

			/* Optical flow on  the entire image is too noise prone
            if (mMOP2fptsPrev.rows () == 0) {
            
                // first time through the loop so we need prev and this mats
                // plus prev points
                // Convert color to gray
                Imgproc.cvtColor (rgbaMat, matOpFlowThis, Imgproc.COLOR_RGBA2GRAY);
                                
                // copy that to prev mat
                matOpFlowThis.copyTo (matOpFlowPrev);
                                
                // get prev corners
                Imgproc.goodFeaturesToTrack (matOpFlowPrev, MOPcorners, iGFFTMax, 0.05, 20);
                mMOP2fptsPrev.fromArray (MOPcorners.toArray ());
                                
                // get safe copy of this corners
                mMOP2fptsPrev.copyTo (mMOP2fptsSafe);
            } else {
                // we've been through before so
				// this mat (matOpFlowThis) is valid. Copy it to prev mat
                matOpFlowThis.copyTo (matOpFlowPrev);
                                
                // get this mat
                Imgproc.cvtColor (rgbaMat, matOpFlowThis, Imgproc.COLOR_RGBA2GRAY);
                                
                // get the corners for this mat
                Imgproc.goodFeaturesToTrack (matOpFlowThis, MOPcorners, iGFFTMax, 0.05, 20);
                mMOP2fptsThis.fromArray (MOPcorners.toArray ());
                                
                // retrieve the corners from the prev mat
                // (saves calculating them again)
                mMOP2fptsSafe.copyTo (mMOP2fptsPrev);
                                
                // and save this corners for next time through
                                
                mMOP2fptsThis.copyTo (mMOP2fptsSafe);
            }

            Video.calcOpticalFlowPyrLK (matOpFlowPrev, matOpFlowThis, mMOP2fptsPrev, mMOP2fptsThis, mMOBStatus, mMOFerr);
            
            if (!mMOBStatus.empty ()) { //overlay the optical flow vectors
                List<Point> cornersPrev = mMOP2fptsPrev.toList ();
                List<Point> cornersThis = mMOP2fptsThis.toList ();
                List<byte> byteStatus = mMOBStatus.toList ();
            
                int x = 0;
                int y = byteStatus.Count - 1;
                                                
                for (x = 0; x < y; x++) {
                    if (byteStatus [x] == 1) {
                        Point pt = cornersThis [x];
                        Point pt2 = cornersPrev [x];
                                    
                        Imgproc.circle (rgbaMat, pt, 5, colorRed, iLineThickness - 1);
                                    
                        Imgproc.line (rgbaMat, pt, pt2, colorRed, iLineThickness);
                    }
                }
            }
            
            //Imgproc.putText (rgbaMat, "W:" + rgbaMat.width () + " H:" + rgbaMat.height () + " SO:" + Screen.orientation, new Point (5, rgbaMat.rows () - 10), Core.FONT_HERSHEY_SIMPLEX, 1.0, new Scalar (255, 255, 255, 255), 2, Imgproc.LINE_AA, false);
            Utils.matToTexture2D (rgbaMat, texture, webCamTextureToMatHelper.GetBufferColors());
            */

			if (!playerCtrl.Started)
				return;

			Quaternion q = Quaternion.identity;
			playerCtrl.VisualUpdate (new Vector3(), q);
        }
    
        /// <summary>
        /// Raises the disable event.
        /// </summary>
        void OnDisable ()
        {
            webCamTextureToMatHelper.Dispose ();
        }
    }
}