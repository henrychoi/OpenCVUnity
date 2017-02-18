using UnityEngine;
using System.Collections;
using System.Runtime.InteropServices;

public static class iOSTorch 
{
	#region invoking c code
	[DllImport ("__Internal")]
	private static extern void _Init();
	[DllImport ("__Internal")]
	private static extern void _Cleanup();
	
	[DllImport ("__Internal")]
	private static extern void _SetTorch(bool state);
	
	[DllImport ("__Internal")]
	private static extern void _TorchOnWithLevel(float level);
	
	[DllImport ("__Internal")]
	private static extern bool _GetTorch();
	
	[DllImport ("__Internal")]
	private static extern float _GetTorchLevel();
	#endregion
	
	#region constructor
	private static int iOSVersion = 10; // this gets set in the static constructor
	
	/// <summary>
	/// Gets the iOS version so we know what stuff we can run.
	/// Most functions require iOS 4.
	/// Torch intensity requires iOS 6.
	/// </summary>
	static iOSTorch()
	{
//		string[] parts = SystemInfo.operatingSystem.Split(' ');
//		int version = 0;
//		foreach (string part in parts)
//		{
//			int.TryParse(part[0].ToString(), out version);
//			if (version > 0) {
//				break;
//			}
//		}
//		iOSVersion = version;
		if (Application.platform == RuntimePlatform.OSXEditor
			|| Application.platform == RuntimePlatform.WindowsEditor) {
			iOSVersion = 0;
			//Debug.Log ("Using Unity Remote");
		}
		//Debug.Log ("iOSVersion = " + iOSVersion);
	}
	#endregion
	
	#region public functions
	/// <summary>
	/// Sets the specified state of the torch.
	/// </summary>
	/// <param name='state'>the state of the torch. True for on. False for off.</param>
	public static void Set(bool state)
	{
		if (iOSVersion >= 4)
			_SetTorch(state);
	}
	
	/// <summary>
	/// Sets the specified level of the torch's brightness.
	/// </summary>
	/// <param name='level'>Range from 0-1. 0 will turn the torch off. 1 is full brightness.</param>
	/// <remarks>iOS 6 only.</remarks>
	public static void Set(float level)
	{
		if (level <= 0)
			Off();
		else
			On(level);
	}
	
	/// <summary>
	/// Gets the current state of the torch.
	/// </summary>
	/// <returns>bool - true for on, false for off.</returns>
	public static bool Get()
	{
		if (iOSVersion >= 4)
			return _GetTorch();
		return false;
	}
	
	/// <summary>
	/// Gets the current brightness level of the torch
	/// </summary>
	/// <returns>float - range from 0-1.</returns>
	/// <remarks>Requires iOS 6</remarks>
	public static float GetLevel()
	{
		if (iOSVersion >= 6)
			return _GetTorchLevel();
		return 0;
	}
	
	/// <summary>
	/// Turns the torch on with the specified brightness.
	/// </summary>
	/// <param name='level'>(optional. defaults to 1) iOS 6+ only. Range from 0-1. 0 will turn the torch to full dimness. 1 is full brightness.</param>
	/// <remarks>At least for the iphone 4 (iOS 6), full dimness is not very dim. Definitely not off.</remarks>
	public static void On(float level = 1)
	{
		if (iOSVersion >= 6) {
			if (level <= 0) {
				Debug.Log ("iOSTorch.cs: Level (" + level + ") is too low. Changing to min (0).");
				level = 0.001f;
			} else if (level > 1) {
				Debug.Log ("iOSTorch.cs: Level (" + level + ") is too high. Changing to max (1).");
				level = 1;
			}
			_TorchOnWithLevel (level);
		} else if (iOSVersion >= 4) {
			//Debug.Log ("Turning on torch");
			Set (true);
		}
	}
	
	/// <summary>
	/// Turns the torch off
	/// </summary>
	public static void Off()
	{
		//Debug.Log ("Turning off torch");
		Set(false);
	}
	
	/// <summary>
	/// Sets up a capture session so there is no delay when turning the torch on. Not required, but recommended.
	/// </summary>
	public static void Init()
	{
		//Debug.Log ("iOSVersion " + iOSVersion);
		if (iOSVersion >= 4) {
			//Debug.Log ("Initing iOSTorch");
			_Init ();
		}
	}
	
	/// <summary>
	/// Disposes of the capture session setup in Init.
	/// </summary>
	/// <remarks>You can turn the torch on after calling this, but there will be a slight delay.</remarks>
	public static void Cleanup()
	{
		if (iOSVersion >= 4)
			_Cleanup();
	}
	#endregion
	
	#region properties
	/// <summary>
	/// Gets or sets the current state of the torch.
	/// </summary>
	/// <value>
	/// The current state of the torch. True means on, False means off.
	/// </value>
	public static bool CurrentState
	{
		get { return Get(); }
		set { Set(value); }
	}
	
	/// <summary>
	/// Gets or sets the current brightness level of the torch.
	/// </summary>
	/// <value>Set is iOS 6+ only. Range from 0-1. 0 will turn the torch to full dimness. 1 is full brightness.</value>
	/// <remarks>At least for the iphone 4 (iOS 6), full dimness is not very dim. Definitely not off.</remarks>
	public static float Level
	{
		get { return GetLevel(); }
		set { On(value); }
	}
	#endregion
}
