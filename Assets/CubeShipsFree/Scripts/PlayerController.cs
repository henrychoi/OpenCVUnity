using UnityEngine;
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

    public class PlayerController : MonoBehaviour
    {
		enum State {
			Waiting4Still,
			Initializing, //Collecting INS input
			ConfirmingBias,//How big is the residual bias after initialization?
			Updating // normal game playing state
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
		const int F_fixed_update = 50;//default is 50
		const float Ts = 1.0f / F_fixed_update;

		#region sensor stats
		Vector3 w_bias_inb_lb = new Vector3(); //gyro bias

		// accel specific gravity [g] window;
		//Specific acc INCLUDES gravity (i.e. REACTION to gravity NOT subtracted)
		Vector3Window i_inmu_lmu_ens = new Vector3Window (F_fixed_update, null);//"acc"); // [g]
		// gyro [dps] window
		Vector3Window w_inb_lb_ens = new Vector3Window (F_fixed_update, null);// [rad/s]
		//Vector3Window m_inb_lb_ens = new Vector3Window (F_fixed_update, null);// [uT]
		#endregion

		bool started = false;
		public bool Started { get { return started; } }

		public Text statusText;


		//The current position and rotation of the player in local tangtent frame (l)
		Vector3 r_inl_lb = new Vector3(),
			r_inU_lb  = new Vector3(), //same vector but scaled and rotated from physics to Unity
			g_inl = Vector3.down,
			f_inb_lb = new Vector3(), //includes reaction to gravity, in body frame
			a_inl_lb, //reaction to gravity taken out
			v_inl_lb;

		Vector3Window r_inU_lb_ens = new Vector3Window (F_fixed_update, null); 

		static Quaternion Q_lb = new Quaternion (0, 0, Mathf.Sin (Mathf.PI / 4), Mathf.Cos (Mathf.PI / 4))
		                         * new Quaternion (Mathf.Sin (Mathf.PI / 4), 0, 0, Mathf.Cos (Mathf.PI / 4))
			//, Q2_lb = new Quaternion (0, Mathf.Sin (-Mathf.PI / 4), 0, Mathf.Cos (-Mathf.PI / 4))
			//	* new Quaternion (Mathf.Sin (Mathf.PI / 2), 0, 0, Mathf.Cos (Mathf.PI / 2))
			;
		Quaternion q_lb = Q_lb, q_lb_avg = Q_lb //, q_bl_avg = Quaternion.Inverse(Q_lb)
			, q_Ub = Quaternion.identity
			, q_inc_lb = Quaternion.identity
		;

		static Quaternion Q_Ul = //Will be used often, so calculate once
			new Quaternion(0, Mathf.Sin(-Mathf.PI/4), 0, Mathf.Cos(-Mathf.PI/4))
			* new Quaternion(Mathf.Sin(-Mathf.PI/4), 0, 0, Mathf.Cos(-Mathf.PI/4))
			;

		public void VisualUpdate(Vector3 dr, Quaternion dtheta) {
			//Cutting the tie from visual odometry (not that optical flow is visual odometry)
//			r += dr * Scale_w2game;
//			q *= dtheta;
//			Rternion.Normalize (ref q);
		}

		// Need a handle to the PiP camera transform, to unrotate it back to inertial frame
		public Camera attitudeCam, positionCam;//These are connected to this class in the inspector
		Vector3 attitudeCamOffset, positionCamOffset;
		bool yIsUp = true;//alternative is up

		// Use this for initialization
        void Start()
		{
			#if WHAT_IS_UNITY_QUATERNION_MULTIPLY_DOING //Understand Unity's quaternion rotation
			Quaternion q = Q_lb
				//new Quaternion (0, Mathf.Sin (-Mathf.PI / 4), 0, Mathf.Cos (-Mathf.PI / 4)) *
			    //new Quaternion (Mathf.Sin (-Mathf.PI / 4), 0, 0, Mathf.Cos (-Mathf.PI / 4))
				;
			Vector3 v = q * Vector3.up;
			Matrix4x4 mx = Matrix4x4.TRS (Vector3.zero, q, Vector3.one);
			Debug.Log(String.Format("quat {0} mx {2}, rot {1}", q, v, mx));
			//Assert.AreEqual(v, Vector3.forward);
			Assert.IsTrue(false);
			#endif
			//Debug.Log ("platform =" + Application.platform);
			if (Application.platform == RuntimePlatform.OSXEditor
			    || Application.platform == RuntimePlatform.WindowsEditor) {
				yIsUp = true;
				Debug.Log ("Using Unity Remote");
			}

			Input.compensateSensors = true;
			//Must be faster or equal rate than FixedUpdate interval
			Input.gyro.enabled = true;
			//Input.compass.enabled = true;
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
	    //const float HALF_PI = 0.5f * Mathf.PI;
	    //int case_num = 0;

		public void Reset()
		{
			//statusText.text = "Please hold still";
			statusText.gameObject.SetActive(true);
			nZVU = 0;
			state = State.Initializing;
			//Input.compass.enabled = true;
		}

        void FixedUpdate()
        {   //Debug.Log("dT = " + Time.deltaTime);
			Vector3 i_inmu_lmu = Input.acceleration; //sensor measures INERTIAL force 
			Vector3 w_inb_lb = //Input.gyro.rotationRate // in mu frame (RH!)
					Input.gyro.rotationRateUnbiased
				;
			
			//Take from mu to b, going from RH to LF system
			if (yIsUp) {
				w_inb_lb.Set (-w_inb_lb.x, -w_inb_lb.y, w_inb_lb.z);//R-> rotation, so sign change
			} else {
				w_inb_lb.Set (w_inb_lb.y, -w_inb_lb.x, w_inb_lb.z);
			}

			i_inmu_lmu_ens.Add (i_inmu_lmu); //These 2 stats are always useful for initializing
			w_inb_lb_ens.Add (w_inb_lb);     //or updating
			//m_ens.Add(Input.compass.rawVector);
			//float heading = Input.compass.magneticHeading;//[deg]

			switch (state) {
			case State.Waiting4Still:
				if (++nZVU < F_fixed_update)
					break;
				nZVU = 0;
				state = State.Initializing;
				break;
			case State.Initializing:
				//m_inb_lb_ens.Add(Input.compass.rawVector);//just prime the averaging queue
				if (++nZVU < F_fixed_update)
					break;
				//Debug.Log ("heading " + Input.compass.magneticHeading);
				Debug.Log ("i_inmu_lmu " + i_inmu_lmu_ens.Average);
				//Input.compass.enabled = false; //turn off compass to save some power
				w_bias_inb_lb = w_inb_lb_ens.Average;
				Debug.Log (String.Format ("Gyro bias ({0}, {1}, {2})",
					w_bias_inb_lb.x, w_bias_inb_lb.y, w_bias_inb_lb.z));

				//Coarse tip/tilt estimate, see Groves 2nd ed. section 5.6, p. 198
				Vector3 g_mu = i_inmu_lmu_ens.Average.normalized; //Zero acceleration assumption

				//Since I don't (want to) know the latitude and longitude, I can't look up the 
				//world gravity model.  Instead, I will use this as the gravity estimate in the
				g_inl.Set (0, 0, -i_inmu_lmu_ens.Average.magnitude); //Updating state below
				//g_inl.Set(0, 0, -1);

				//Calculate the nominal case
				float roll, pitch;
				if (yIsUp) {
					roll = Mathf.Asin (-g_mu.x);
					pitch = Mathf.Atan2 (-g_mu.z, -g_mu.y);
				} else {
					roll = Mathf.Asin (g_mu.y);
					pitch = Mathf.Atan2 (-g_mu.z, -g_mu.x);
				}
				//statusText.text = String.Format ("roll {0:N1} pitch {1:N1}, a_b {2}"
				//	, roll, pitch, a_b);

				//Assert.IsFalse(case_num == 2);
				//case_num = 1;

				//Reset quaternion to yaw degrees around the z axis,
				//I don't use Quaternion.Euler(Mathf.Rad2Deg * roll, Mathf.Rad2Deg * pitch, 0)
				//to avoid, thinking about left-handed rotation.
				roll *= 0.5f; // Halve the Euler angles feed to quaternion construction
				pitch *= 0.5f;

				Quaternion q = // Apply the roll first, then pitch
					new Quaternion (0, 0, Mathf.Sin (roll), Mathf.Cos (roll))
					* new Quaternion (Mathf.Sin (pitch), 0, 0, Mathf.Cos (pitch));
				//q_lb = Q_lb
				//	* Quaternion.Euler(pitch, 0, roll);//minor miracle: the order matches what I need
				//Debug.Log (String.Format ("roll {0}, pitch {1} a {2} q {3}", roll, pitch, a_b, q));

				const float THRESHOLD = 0.8f;
				if (Mathf.Abs (g_mu.x) < THRESHOLD) { // 1st case: cos(roll) reasonable value
					q_lb = Q_lb * q;
					//statusText.text = String.Format ("roll {0:N1} pitch {1:N1}, a_b {2}", roll, pitch, g_mu);
				} else { // 2nd case, roll near +/- 90; ay << ax
					roll = yIsUp
						? -Mathf.Sign (g_mu.x) * Mathf.Acos (-g_mu.y)
						: Mathf.Sign (g_mu.y) * Mathf.Acos (-g_mu.x);
					float yaw = Mathf.Asin (g_mu.z / Mathf.Sin (roll));
					//statusText.text = String.Format ("roll {0:N1} yaw {1:N1}, a_b {2}", roll, yaw, g_mu);

					//Assert.IsFalse(case_num == 1);
					//case_num = 2;

					//Reset quaternion to yaw degrees around the z axis,
					//I don't use Quaternion.Euler(Mathf.Rad2Deg * roll, Mathf.Rad2Deg * pitch, 0)
					//to avoid, thinking about left-handed rotation.
					roll *= 0.5f; // Halve the Euler angles feed to quaternion construction
					yaw *= 0.5f;

					Quaternion q2 = // Apply the roll first, then yaw about Y
						new Quaternion (0, 0, Mathf.Sin (roll), Mathf.Cos (roll))
						* new Quaternion (0, Mathf.Sin (yaw), 0, Mathf.Cos (yaw));
					//Debug.Log (String.Format("roll {0}, yaw {1} a {2} q {3}", roll, yaw, a_b, q_lb));

					q_lb = Q_lb * Quaternion.Slerp (q, q2, Mathf.Abs (g_mu.x) - THRESHOLD);
				}
				//Sanity check the solution: q_ub(i_inmu_lmu_ens.Average) rotated by q_bl should be
				//g_inl, but since Quaternion.operator*(Quaternion q, Vector3 v) actually rotates
				//the vector the OTHER way, do some mental juggling
				Vector3 i_b_avg = i_inmu_lmu_ens.Average;
				i_b_avg.Set (i_b_avg.x, i_b_avg.y, -i_b_avg.z); // IMU --> b frame
				//Debug.Log ("g_inmu " + g_inmu);
				Vector3 i_inl = q_lb * i_b_avg;
				Assert.AreApproximatelyEqual (i_inl.x, g_inl.x);
				Assert.AreApproximatelyEqual (i_inl.y, g_inl.y);
				Assert.AreApproximatelyEqual (i_inl.z, g_inl.z);

				//Coarse heading estimate, according to Groves 2nd ed. section 6.1.1, p. 219
//				float s_roll = Mathf.Sin (roll), c_roll = Mathf.Cos (roll)
//				, s_pitch = Mathf.Sin (pitch), c_pitch = Mathf.Cos (pitch)
//			    , m_num = -m_b.y * c_roll + m_b.z * s_roll
//				, m_den = m_b.x + m_b.y + m_b.z;

				r_inl_lb.Set (0, 0, 0);//Q: shall I reset the position at 1 km height?
				v_inl_lb.Set (0, 0, 0);
				//statusText.gameObject.SetActive(false);
				state = State.Updating;//State.ConfirmingBias;
				nZVU = 0;
				break;

			case State.ConfirmingBias:
				if (++nZVU < F_fixed_update)
					break;
				Vector3 w_inb = w_inb_lb_ens.Average;
			    i_b_avg = i_inmu_lmu_ens.Average;
				Assert.AreApproximatelyEqual(i_b_avg.x, i_inmu_lmu_ens.Average.x);
				Assert.AreApproximatelyEqual(i_b_avg.y, i_inmu_lmu_ens.Average.y);
				Assert.AreApproximatelyEqual(i_b_avg.z, i_inmu_lmu_ens.Average.z);

				state = State.Updating;
				break;

			default: //navigation update
				//Debug.Log ("a_inmu_lmu " + a_inmu_lmu);
				//Assert.IsTrue (a_inmu_lmu.magnitude < 2);
				//Debug.Log(String.Format("a {0}, w {1}, gyro {2}", a_inb_lb, w_inb_lb, Input.gyro.enabled));
				//q_inc_lb.Set (0, 0, 0.5f * 0.01f, 1);//For now, let's make up a rotation vector

				//Calculate the average attitude (to rotate the specific force to l)
				//See Groves GNSS 2nd ed. Appendix E section 6.3
				Vector3 alpha = 0.5f * Time.deltaTime * w_inb_lb;
				float As = 0.5f, Ac = 1.0f; //1st order approximation
				const int O_ATTITUDE_UPDATE = 4;
				float alpha_div2, alpha_div2_sq, alpha_div2_qd;
				if (O_ATTITUDE_UPDATE == 4) {
					alpha_div2 = 0.5f * alpha.magnitude;
					alpha_div2_sq = alpha_div2 * alpha_div2;
					alpha_div2_qd = alpha_div2_sq * alpha_div2_sq;
					As = 0.5f - 0.083333333333333f * alpha_div2_sq;
					Ac = 1.0f - 0.5f * alpha_div2_sq + 0.041666666666667f * alpha_div2_qd;
				}
				alpha *= As;
				q_inc_lb.Set (alpha.x, alpha.y, alpha.z, Ac); // 4st order approximation
				q_lb_avg = q_lb * q_inc_lb;//Rotate the current attitude by (attitude increment)/2

				//No need to invert the q_lb to q_bl, because q_lb * v actually works like rotate(q_bl, v)
				//q_bl_avg = Quaternion.Inverse (q_lb_avg);
				//q_bl_avg.Set (q_lb_avg.x, q_lb_avg.y, q_lb_avg.z, -q_lb_avg.w);

				//Specific force on body is equal and opposite to INERTIAL reaction,
				//but Z flips sign when going R->L frame //Still has REACTION to gravity
				f_inb_lb.Set (-i_inmu_lmu.x, -i_inmu_lmu.y, i_inmu_lmu.z);
				//Rotate force in b to l and ADD gravity (because f_inb had REACTION to gravity)
				const float G2MPS2 = 9.807f; //from wiki
				//Vector3 f_inl_lb = q_bl_avg * f_inb_lb;
				a_inl_lb = G2MPS2//* (f_inl_lb + g_inl);//
					* (q_lb_avg * f_inb_lb + g_inl);//* (q_bl_avg * f_inb_lb + g_inl);
				//a_inl_lb = new Quaternion (.5f, .5f, .5f, -.5f) * Vector3.up;
				Vector3 v_inl_lb_p = v_inl_lb + a_inl_lb * Time.deltaTime;
				Vector3 v_inl_avg = 0.5f * (v_inl_lb + v_inl_lb_p);
				v_inl_lb = v_inl_lb_p; //update to new velocity
				r_inl_lb += v_inl_avg * Time.deltaTime;
					
				statusText.text = String.Format ("f_inb {2} a_inl {0}, v_inl_avg {1}",
					a_inl_lb, v_inl_avg, f_inb_lb);
				//Debug.Log(String.Format("a_inl_lb {0}, v_inl_lb_avg {1}", a_inl_lb, v_inl_avg));
				//Assert.IsTrue (a_inl_lb.magnitude < 2.0f);//sanity test

				//Compute the new attitude (not the average)
				alpha = Time.deltaTime * w_inb_lb;
				alpha_div2 = 0.5f * alpha.magnitude;
				alpha_div2_sq = alpha_div2 * alpha_div2;
				alpha_div2_qd = alpha_div2_sq * alpha_div2_sq;
				As = 0.5f - 0.083333333333333f * alpha_div2_sq;
				Ac = 1.0f - 0.5f * alpha_div2_sq + 0.041666666666667f * alpha_div2_qd;
				alpha *= As;
				q_inc_lb.Set (alpha.x, alpha.y, alpha.z, Ac); // 4st order approximation
				//q_inc_lb.Set (0, 0, 0.5f * 0.01f, 1);//make up a rotation for debugging
				q_lb = q_lb * q_inc_lb;//Rotate the current attitude by the attitude increment

				break;
			}

			// Update the Unity position and attitude
			r_inU_lb.Set(-r_inl_lb.y, r_inl_lb.z, r_inl_lb.x);
			const int Scale_l2U = 1000/10;// physics is in 10 cm scale (phone dimension); Unity is in m
			r_inU_lb *= Scale_l2U;

			q_Ub = Q_Ul * q_lb; //Rotate the q_lb by Q_lU to go to Unity
			//q_Ub.Set(-q_Ub.x, -q_Ub.z, -q_Ub.y, q_Ub.w); //then go to Unity frame (LH)
			Rternion.Normalize(ref q_Ub);
			//Debug.Log("Moving to " + this.r);
			myRigidbody.MovePosition(r_inU_lb); myRigidbody.MoveRotation(q_Ub);
			//myRigidbody.position = r; myRigidbody.rotation = q;

			r_inU_lb_ens.Add(r_inU_lb);

			//Move the chase cameras (shown in PiP) with the player
			attitudeCam.transform.position = r_inU_lb + attitudeCamOffset;
			positionCam.transform.position = r_inU_lb_ens.Average + positionCamOffset;
        }

        void Update()
        {
			RaycastHit hitInfo;
			if (Physics.Raycast(r_inl_lb //The starting point of the ray in world coordinates.
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
