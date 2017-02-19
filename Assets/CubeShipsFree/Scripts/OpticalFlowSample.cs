using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;
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
    
        // Use this for initialization
        void Start ()
        {
            webCamTextureToMatHelper = gameObject.GetComponent<WebCamTextureToMatHelper> ();
            webCamTextureToMatHelper.Init ();
    
			//gameObject is what this script is attached to.  In this project, it is the Quad
			//Quad's parent is the Player, which has the PlayerController script component
			playerCtrl = gameObject.GetComponentInParent<PlayerController> ();
			Assert.IsNotNull (playerCtrl);
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

			if (rxgxbMat != null) rxgxbMat.Dispose ();
			if (tempA != null) tempA.Dispose ();
			if (tempB != null) tempB.Dispose ();

			if (matOpFlowThis != null)
                matOpFlowThis.Dispose ();
            if (matOpFlowPrev != null)
                matOpFlowPrev.Dispose ();
            if (MOPcorners != null)
                MOPcorners.Dispose ();
            if (mMOP2fptsThis != null)
                mMOP2fptsThis.Dispose ();
            if (mMOP2fptsPrev != null)
                mMOP2fptsPrev.Dispose ();
            if (mMOP2fptsSafe != null)
                mMOP2fptsSafe.Dispose ();
            if (mMOBStatus != null)
                mMOBStatus.Dispose ();
            if (mMOFerr != null)
                mMOFerr.Dispose ();
        }

        /// <summary>
        /// Raises the web cam texture to mat helper error occurred event.
        /// </summary>
        /// <param name="errorCode">Error code.</param>
        public void OnWebCamTextureToMatHelperErrorOccurred(WebCamTextureToMatHelper.ErrorCode errorCode){
            Debug.Log ("OnWebCamTextureToMatHelperErrorOccurred " + errorCode);
        }

		float lastOdoTime = .0f;

		Mat rxgxbMat, tempA, tempB,
			emptyMat = new Mat(),
			erode_kernel = new Mat()
				//Imgproc.getStructuringElement (Imgproc.MORPH_ELLIPSE, new Size(7,7))
		;
		Point erode_anchor = new Point(0,0);

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

			// Assume the brightest white region(s) are the torches
			List<Mat> channels = new List<Mat>();
			Core.split (rgbaMat, channels);
			//Debug.Log ("channels: " + channels.Count);
			if (tempA == null) tempA = new Mat (height, width, channels [0].type ());
			if (tempB == null) tempB = new Mat (height, width, channels [0].type ());

			//Core.mulSpectrums(channels[0], channels[1], rxgxbMat

			rxgxbMat = channels[0].mul(channels[1].mul(channels[2], 1.0f/256), 1.0f/256);
			Imgproc.medianBlur (rxgxbMat, tempA, 7);
			//Imgproc.erode (tempA, tempB, erode_kernel, erode_anchor, 2);
			Imgproc.morphologyEx (tempA, tempB, Imgproc.MORPH_OPEN, erode_kernel, erode_anchor, 2);

			/*
			Imgproc.calcHist (new List<Mat> (new Mat[]{ rxgxbMat }),
			    new MatOfInt (0), emptyMat,empty maskMat => no mask
				histMat, //output
				new MatOfInt (16), //histSize(s): # bins
				new MatOfFloat (0, 180)//the left and right most bins
			);
			Debug.Log ("hist " + histMat);
			*/

			//Debug.Log("max of cross image = " + rxgxbMat.
			//Debug.Log (String.Format ("original size {0}x{1}", rgbaMat.width(), rgbaMat.height()));
			//Debug.Log (String.Format ("new size {0}x{1}", rxgxbMat.width(), rxgxbMat.height()));

			Utils.matToTexture2D (tempB, texture, webCamTextureToMatHelper.GetBufferColors());

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