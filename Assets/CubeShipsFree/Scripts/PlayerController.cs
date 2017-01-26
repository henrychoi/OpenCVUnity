using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System;

// Heavily based on Unity Space Shooter tutorial
namespace CubeSpaceFree
{

    [System.Serializable]
    public class Boundary
    {
        public float xMin, xMax, zMin, zMax;
    }

	public class Vector3Window {
		Queue<Vector3> q1, q2;
		Vector3 average = new Vector3(), variance = new Vector3(), stddev = new Vector3();
		public Vector3 Average { get { return average; } }
		public Vector3 Variance { get { return variance; } }
		public Vector3 StdDev { get { return stddev; } }
		int size;
		public int Size { get { return size; } }
		float sizeinv;
		public bool Valid { get { return q1.Count == size; } }

		public Vector3Window(int sz) {
			q1 = new Queue<Vector3> (sz);
			q2 = new Queue<Vector3> (sz);
			size = sz;
			sizeinv = 1.0f / size;
		}
		public void Add(Vector3 v) {
			if (q1.Count >= size) { // full, so have to evict the earliest
				Vector3 head = q1.Dequeue();
				average -= head;
				head = q2.Dequeue ();
				variance -= head;
			}

			Vector3 scaled = v * sizeinv;
			q1.Enqueue (scaled);
			average += scaled; // simple accumulation

			// Running variance += (x[[k]/size) * x[k]
			scaled.x *= v.x; scaled.y *= v.y; scaled.z *= v.z;
			variance += scaled;
			q2.Enqueue (scaled);

			stddev.x = Mathf.Sqrt (variance.x);
			stddev.y = Mathf.Sqrt (variance.y);
			stddev.z = Mathf.Sqrt (variance.z);
		}
	}

    public class PlayerController : MonoBehaviour
    {
		#region MOVEMENT CONTROL
        //public float speed, tilt;              // tilt factor
        //public Boundary boundary;       // movememnt boundary
		//public float smoothing = 5;     // this value is used for smoothing movement
		//private Vector3 smoothDirection;// used to smooth out mouse and touch control
		public Rigidbody myRigidbody;//Necessary for physics; will the player be subject to physics?
		#endregion

        public GameObject shot;         // bullet prefab
        public Transform shotSpawn;     // the turret (bullet spawn location)
        public float fireRate = 0.5f;   
        private float nextFire = 0.0f;

		#region sensor stats
		//Vector3 a = new Vector3(), w = new Vector3();
		Vector3Window a_ens = new Vector3Window (25); // 0.5 sec accel window
		Vector3Window w_ens = new Vector3Window (25); // 0.5 sec gyro window
		#endregion

		bool started = false;
		public bool Started { get { return started; } }

		//The current position and rotation of the player
		Vector3 r = new Vector3();//Unity’s default unit scale is 1 unit = 1 meter
		Quaternion q = new Quaternion {x=0, y=0, z=0, w=1};
		Vector3Window r_ens = new Vector3Window (25); // 0.5 sec accel window// low passed position

		void NormalizeRotQ(ref Quaternion q) {
			float l2inv = 1.0f/Mathf.Sqrt (q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
			q.x *= l2inv;
			q.y *= l2inv;
			q.z *= l2inv;
			q.w *= l2inv;
		}

		public void VisualUpdate(Vector3 dr, Quaternion dtheta) {
			//Debug.Log ("visual update dr = " + dr);
			r += dr;
			q *= dtheta;
			NormalizeRotQ (ref q);
		}

		// Need a handle to the PiP camera transform, to unrotate it back to inertial frame
		public Camera attitudeCam, positionCam;//These are connected to this class in the inspector
		Vector3 attitudeCamOffset, positionCamOffset;

        // Use this for initialization
        void Start()
		{
			Input.compensateSensors = false;
			//Would like 100 Hz, but accel is ticking at 50 Hz anyway; consider increasing later
			//if I need to update the attitude more frequently than integrating the acceleration.
			Input.gyro.updateInterval = 0.02f;

			attitudeCamOffset = attitudeCam.transform.position;
			positionCamOffset = positionCam.transform.position;

            myRigidbody = GetComponent<Rigidbody>();

			started = true;
        }


		//Vector3 origin = new Vector3();
		//void OnDrawGizmos() { //This IS called, but I just did not see the ray
		//	Debug.Log ("I am here");
		//	Debug.DrawRay(r, a, Color.yellow, 1.0f, false);
			//Ray ray = new Ray {}
			//Gizmos.DrawLine(origin, new Vector3 {x=1, y=0, z=0});
		//}

        void FixedUpdate()
        {
			//Debug.Log("dT = " + Time.deltaTime);

			#if KEYBOARD_MOUSE_INPUT
            // keyboard
            float moveHorizontal = Input.GetAxis("Horizontal");
            float moveVertical = Input.GetAxis("Vertical");

            // if keyboard direction key is pressed
            if (moveHorizontal != 0 || moveVertical != 0)
            {
                myRigidbody.velocity = new Vector3(moveHorizontal, 0.0f, moveVertical) * speed;

            }
            else
            {
                Vector3 pos = Input.mousePosition;

                pos.z = Camera.main.transform.position.y + 1;
                pos = Camera.main.ScreenToWorldPoint(pos);
                Vector3 origin = new Vector3(transform.position.x, transform.position.y, transform.position.z);//.zero;

                Vector2 currentPosition = new Vector3(pos.x, pos.z);
                Vector3 directionRaw = pos - origin;
                //Debug.Log("directionRaw.magnitude=" + directionRaw.magnitude);

                Vector3 direction = directionRaw.normalized;

                smoothDirection = Vector3.MoveTowards(smoothDirection, direction, smoothing);

                direction = smoothDirection;
                Vector3 movement = new Vector3(direction.x, 0, direction.z);
                myRigidbody.velocity = movement * speed;
                Debug.Log("currentPosition=" + currentPosition + "  Input.mousePosition=" + Input.mousePosition +
                    " pos =" + pos + " direction=" + direction + " smoothDirection=" + smoothDirection +
                    " movement=" + movement);

                //transform.position = Vector3.Lerp(transform.position, pos, Time.deltaTime* speed*2);
                //myRigidbody.velocity = new Vector3(moveHorizontal, 0.0f, moveVertical) * speed;                    
            }

            transform.position = new Vector3
            (
                Mathf.Clamp(transform.position.x, boundary.xMin, boundary.xMax),
                0.0f,
                Mathf.Clamp(transform.position.z, boundary.zMin, boundary.zMax)
            );
            myRigidbody.rotation = Quaternion.Euler(0, 0, myRigidbody.velocity.x * -tilt);

			#else
			//not sure if this is compensated, but still estimate bias
			a_ens.Add(Input.acceleration);
			//Debug.Log(String.Format("acc = ({0:F},{1:F},{2:F})", acc.x, acc.y, acc.z));

			//I will do further bias estimation using the compensated gyro values
			w_ens.Add(Input.gyro.rotationRateUnbiased);

			//Just use 6-axis for now
			//float heading = Input.compass.magneticHeading;//[deg]

			//Debug.Log("Moving to " + this.r);
			myRigidbody.MovePosition(r); myRigidbody.MoveRotation(q);
			//myRigidbody.position = r; myRigidbody.rotation = q;

			r_ens.Add(r);

			//Move the chase cameras (shown in PiP) with the player
			attitudeCam.transform.position = r + attitudeCamOffset;
			positionCam.transform.position = r_ens.Average + positionCamOffset;
			#endif
        }

		/*Don't fire a bullet in this navigation demo app
        void Update()
        {
            if ((Input.GetButton("Fire1") || Input.GetKeyDown(KeyCode.Space)) && Time.time > nextFire)
            {
                nextFire = Time.time + fireRate;
                Instantiate(shot, shotSpawn.position, shotSpawn.rotation);
                GetComponent<AudioSource>().Play();
            }
        }
        */
    }
}
