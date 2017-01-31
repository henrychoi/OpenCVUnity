using System;
using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Rector3 {
	public static float x_l(Vector3 v) { return v.x; }
	public static float y_l(Vector3 v) { return v.y; }
	public static float z_l(Vector3 v) { return -v.z; }
	public static void toLeft(ref Vector3 l, Vector3 v) {
		l.Set(v.x, v.y, -v.z);
	}
}

//@brief To convert a RHS Q to LHS: (-x,-z,-y, w)
public class Rternion  {
	public float x_l(Quaternion q) { return -q.x; } 
	public float y_l(Quaternion q) { return -q.z; } 
	public float z_l(Quaternion q) { return -q.y; } 
	public float w_l(Quaternion q) { return q.w; }
	//public Quaternion left(Quaternion q) {
	//	return new Quaternion(-q.x, -q.z, -q.y, q.w);
	//}
	public void toLeft(ref Quaternion l, Quaternion q) {
		l.Set(-q.x, -q.z, -q.y, q.w);
	}
	public static void Normalize(ref Quaternion q) {
		float l2inv = 1.0f/Mathf.Sqrt (q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
		q.x *= l2inv;
		q.y *= l2inv;
		q.z *= l2inv;
		q.w *= l2inv;
	}
}


public class Vector3Window {
	Queue<Vector3> q1, q2;
	Vector3 m1 = new Vector3()
		, m2 = new Vector3()
		, variance = new Vector3(), stddev = new Vector3();
	public Vector3 Average { get { return m1; } }
	public Vector3 Variance { get { return variance; } }
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