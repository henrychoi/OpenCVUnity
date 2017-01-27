using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine.Assertions;
using UnityEngine.UI;

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
		Vector3 m1 = new Vector3()
			, m2 = new Vector3()
			, variance = new Vector3(), stddev = new Vector3();
		public Vector3 Average { get { return m1; } }
		public Vector3 Variance { get { return m2; } }
		public Vector3 StdDev { get { return stddev; } }
		int size;
		public int Size { get { return size; } }
		float sizeinv;
		public bool Valid { get { return q1.Count == size; } }
		string debugStr;

		public Vector3Window(int sz, string name) {
			q1 = new Queue<Vector3> (sz);
			q2 = new Queue<Vector3> (sz);
			size = sz;
			sizeinv = 1.0f / size;
			this.debugStr = name;
		}
		public void Add(Vector3 v) {
			if (q1.Count >= size) { // full, so have to evict the earliest
				Vector3 head = q1.Dequeue();
				m1 -= head;
				head = q2.Dequeue ();
				m2 -= head;
			}

			Vector3 scaled1 = v * sizeinv;
			q1.Enqueue (scaled1);
			m1 += scaled1; // simple accumulation

			// Running average variance += (x[[k]/size) * x[k] - mean^2
			Vector3 scaled2 = scaled1;
			scaled2.x *= v.x; scaled2.y *= v.y; scaled2.z *= v.z;
			q2.Enqueue (scaled2);
			m2 += scaled2;

			variance.Set(Mathf.Max(m2.x - m1.x * m1.x, 0)
				, Mathf.Max(m2.y - m1.y * m1.y, 0)
				, Mathf.Max(m2.z - m1.z * m1.z, 0));

			stddev.x = Mathf.Sqrt (variance.x);
			stddev.y = Mathf.Sqrt (variance.y);
			stddev.z = Mathf.Sqrt (variance.z);

			if (debugStr != null)
				Debug.Log (String.Format ("{0} {1}/{2} <- {3}: {4}", debugStr, m1, m2, v, variance));
		}
	}

    public class PlayerController : MonoBehaviour
    {
		enum State {
			Resetting,
			Updating
		}
		State state = State.Updating;
		int nZVU;

		#region MOVEMENT CONTROL
        //public float speed, tilt;              // tilt factor
        //public Boundary boundary;       // movememnt boundary
		//public float smoothing = 5;     // this value is used for smoothing movement
		//private Vector3 smoothDirection;// used to smooth out mouse and touch control
		public Rigidbody myRigidbody;//Necessary for physics; will the player be subject to physics?
		#endregion

		public GameObject enemy;
		TextMesh enemy_text_mesh;

        public GameObject shot;         // bullet prefab
        Transform turret_xform;// the turret (bullet spawn location)
        public float fireRate = 0.5f;   
        // float nextFire = 0.0f;

		// Changed in Unity editor Edit -> Project Settings -> Time -> Fixed Timestep
		const int F_fixed_update = 100;//default is 50

		// physics is in world scale; 
		const int Scale_w2game = 100;

		#region sensor stats
		Vector3 w_bias = new Vector3(); //gyro bias
		Vector3Window a_ens//Specific acc INCLUDES gravity (i.e. REACTION to gravity NOT subtracted)
			= new Vector3Window ((int)(1.0f * F_fixed_update), null) // accel [g] window
			, w_ens = new Vector3Window ((int)(1.0f * F_fixed_update), null)// gyro [dps] window
			//, m_ens = new Vector3Window (0.5 * F_fixed_update) // magnetic [uT] field window
			;
		#endregion

		bool started = false;
		public bool Started { get { return started; } }

		public Text statusText;


		//The current position and rotation of the player
		Vector3 r = new Vector3();//Unity’s default unit scale is 1 unit = 1 meter
		Quaternion q = new Quaternion {x=0, y=0, z=0, w=1};
		Vector3Window r_ens = new Vector3Window ((int)(0.5f * F_fixed_update), null); 

		void NormalizeRotQ(ref Quaternion q) {
			float l2inv = 1.0f/Mathf.Sqrt (q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
			q.x *= l2inv;
			q.y *= l2inv;
			q.z *= l2inv;
			q.w *= l2inv;
		}

		public void VisualUpdate(Vector3 dr, Quaternion dtheta) {
			//Debug.Log ("visual update dr = " + dr);
			r += dr * Scale_w2game;
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

			//Must be faster or equal rate than FixedUpdate interval
			Input.gyro.updateInterval = 1.0f/F_fixed_update;

			//Chase camera position offsets
			attitudeCamOffset = attitudeCam.transform.position;
			positionCamOffset = positionCam.transform.position;

            myRigidbody = GetComponent<Rigidbody>();
			Assert.IsNotNull (myRigidbody);

			turret_xform = transform.Find ("Turret");
			Assert.IsNotNull (turret_xform);

			enemy_text_mesh = enemy.GetComponent<TextMesh> ();
			Assert.IsNotNull (enemy_text_mesh);

			statusText.gameObject.SetActive(false);

			started = true;
        }


		//Vector3 origin = new Vector3();
		//void OnDrawGizmos() { //This IS called, but I just did not see the ray
		//	Debug.Log ("I am here");
		//	Debug.DrawRay(r, a, Color.yellow, 1.0f, false);
			//Ray ray = new Ray {}
			//Gizmos.DrawLine(origin, new Vector3 {x=1, y=0, z=0});
		//}
	    const float HALF_PI = 0.5f * Mathf.PI;
		public void Reset()
		{
			//statusText.text = "Please hold still";
			statusText.gameObject.SetActive(true);
			nZVU = 0;
			state = State.Resetting;
		}

        void FixedUpdate()
        {   //Debug.Log("dT = " + Time.deltaTime);

			//not sure if this is compensated, but still estimate bias
			a_ens.Add(Input.acceleration);
			//Debug.Log(String.Format("acc = ({0:F},{1:F},{2:F})", acc.x, acc.y, acc.z));

			//I will do further bias estimation using the compensated gyro values
			w_ens.Add(Input.gyro.rotationRateUnbiased);

			//m_ens.Add(Input.compass.rawVector);
			//float heading = Input.compass.magneticHeading;//[deg]

			switch(state) {
			case State.Resetting:
				if(++nZVU < F_fixed_update) break;

				w_bias = w_ens.Average;
				//Coarse tip/tilt estimate, see Groves 2nd ed. section 5.6, p. 198
				Vector3 a_b = a_ens.Average //specific force includes gravity
					//, m_b = m_ens.Average
					; // magnetic field
				float roll = Mathf.Atan2 (-a_b.y, -a_b.z)
					, f_yz = Mathf.Sqrt (a_b.y * a_b.y + a_b.z * a_b.z)
					, pitch = Mathf.Atan2 (a_b.x, f_yz) //less work than the explicit form below
					//Mathf.Abs(f_b.x) < f_yz
					//? Mathf.Atan(f_b.x / f_yz) //well-conditioned case
					//: HALF_PI - Mathf.Atan(f_yz / f_b.x)
					;
				Debug.Log (String.Format("resetting with roll {0} pitch {1}, w_bias {2}, Pw_b {3}"
					, roll, pitch, w_bias, w_ens.Variance));

				//Coarse heading estimate, according to Groves 2nd ed. section 6.1.1, p. 219
				//			float s_roll = Mathf.Sin (roll), c_roll = Mathf.Cos (roll)
				//				, s_pitch = Mathf.Sin (pitch), c_pitch = Mathf.Cos (pitch)
				//			    , m_num = -m_b.y * c_roll + m_b.z * s_roll
				//				, m_den = m_b.x + m_b.y + m_b.z;

				statusText.gameObject.SetActive(false);
				state = State.Updating;
				break;

			default: //navigation update
				break;
			}

			//Debug.Log("Moving to " + this.r);
			myRigidbody.MovePosition(r); myRigidbody.MoveRotation(q);
			//myRigidbody.position = r; myRigidbody.rotation = q;

			r_ens.Add(r);

			//Move the chase cameras (shown in PiP) with the player
			attitudeCam.transform.position = r + attitudeCamOffset;
			positionCam.transform.position = r_ens.Average + positionCamOffset;
        }

        void Update()
        {
			RaycastHit hitInfo;
			if (Physics.Raycast(r //The starting point of the ray in world coordinates.
				, turret_xform.forward // direction
				, out hitInfo //hitInfo will contain more information about where the collider was hit
				//, float maxDistance = Mathf.Infinity
				//, int layerMask = DefaultRaycastLayers
				//, QueryTriggerInteraction queryTriggerInteraction = QueryTriggerInteraction.UseGlobal
				)
			) { // can hit the enemy
				GameObject hit = hitInfo.transform.gameObject;
				if (hit == enemy) {
					enemy_text_mesh.color = Color.red;
					return;
				}
			}
			// else will not hit the target
			enemy_text_mesh.color = Color.black;

			/*Don't fire a bullet in this navigation demo app
            if ((Input.GetButton("Fire1") || Input.GetKeyDown(KeyCode.Space)) && Time.time > nextFire)
            {
                nextFire = Time.time + fireRate;
                Instantiate(shot, turret_xform.position, turret_xform.rotation);
                GetComponent<AudioSource>().Play();
            }
			        */
        }
    }
}
