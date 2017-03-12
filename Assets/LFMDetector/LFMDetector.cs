using UnityEngine;
using System.Collections;
using System;
using System.Runtime.InteropServices;

public class LFMDetector
//	: MonoBehaviour
{
	#if UNITY_IOS
	const string dll_name = "__Internal";
	#else
	const string dll_name = "ASimplePlugin";
	#endif

	//Lets make our calls from the Plugin
	[DllImport (dll_name)]
	private static extern int PrintANumber();
	
	[DllImport (dll_name)]
	public static extern IntPtr PrintHello();
	
	[DllImport (dll_name)]
	private static extern int AddTwoIntegers(int i1,int i2);

	[DllImport (dll_name)]
	private static extern float AddTwoFloats(float f1,float f2);	
	
//	void Start () {
//		Debug.Log(PrintANumber());
//		Debug.Log(Marshal.PtrToStringAuto (PrintHello()));
//		Debug.Log(AddTwoIntegers(2,2));
//		Debug.Log(AddTwoFloats(2.5F,4F));
//	}
}
