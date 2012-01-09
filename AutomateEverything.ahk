;
; AutoHotkey Version: 1.x
; Language:       English
; Platform:       Win9x/NT
; Author:         A.N.Other <myemail@nowhere.com>
;
; Script Function:
;	Template script (you can customize this template by editing "ShellNew\Template.ahk" in your Windows folder)
;

Loop, %A_WorkingDir%\* , 2 , 0  
{
	;LOAD PROFUSION
	Run C:\Program Files\Compumedics\ProFusion PSG 2\ProFusionPSG.exe
	WinWait, ProFusion
	Sleep 1000

	;LOAD ith SLEEP STUDY
	;Click 430, 450
	SendInput {Tab}
	Sleep 100
	SendInput {Tab}
	Sleep 100
	SendInput {Space}
	Sleep 500
	SendInput {Tab}
	Sleep 100
	SendInput {Tab}
	Sleep 100
	SendInput {Tab}
	SendInput %A_LoopFileFullPath%
	Sleep 2000
	SendInput {Enter}
	Sleep 3000
	Click 269, 92
	Sleep 500
	SendInput {Enter}

	Sleep 3000

	;EXPORT PATIENT EVENT DATA
	IfWinExist ProFusion
	{ 
		WinActivate
	}

	;Opens Up the Event Window
	SendInput {F12}

	Sleep 600
	;Selects RespEvents
	Click 287,75
	Sleep 600
	SendInput R
	Sleep 600
	IfWinExist Scored
	{
		WinActivate
	}
	Click 677,94
	SendInput {del}{del}{del}
	Sleep 200
	SendInput %A_LoopFileFullPath%\RespEvents_%A_LoopFileName%.csv
	SendInput !S
	Sleep 500
	
	;Selects Arousals
	Click 287,75
	Sleep 500
	SendInput A
	Sleep 500
	IfWinExist Scored
	{
		WinActivate
	}
	Click 677,94
	SendInput {del}{del}{del}
	Sleep 200
	SendInput %A_LoopFileFullPath%\Arousals_%A_LoopFileName%.csv
	Sleep 200
	SendInput !S
	Sleep 500

	;Selects SpO2 Events
	Click 287,75
	Sleep 500
	SendInput S
	Sleep 500
	IfWinExist Scored
	{
		WinActivate
	}
	Click 677,94
	SendInput {del}{del}{del}
	Sleep 200
	SendInput %A_LoopFileFullPath%\SpO2_%A_LoopFileName%.csv
	SendInput !S
	Sleep 200
	Click 673,202
	Sleep 500

	Sleep 1000
	;EXPORT HYPNOGRAM AS TEXT FILE
	IfWinExist ProFusion
	{ 
		WinActivate
	}
	SendInput !E
	Sleep 200
	SendInput H
	Sleep 200
	SendInput E
	Sleep 200
	SendInput %A_LoopFileFullPath%\RK_%A_LoopFileName%.txt
	SendInput {Enter}
	Sleep 200


	;EXPORT PATIENT SIGNALS TO MATLAB
	;Activate the Matlab Add on
	SendInput !O
	Sleep 100
	SendInput M
	Sleep 100
	SendInput {Space}
	SendInput {Tab}{Tab}
	clipboard = ;empty the clipboard
	SendInput ^c
	ClipWait
	clipboard -= 1
	SendInput %clipboard%
	Sleep 200
	Click 50,160
	Sleep 100
	Click 450,60
	;wait for export box to pop up and then press enter
	; pop up window is MatLabExport2
	;Click 125,110 ;this should click ok in new pop up
	WinWait, MatLabExport2
	SendInput {Enter}

	IfWinExist MATLAB
	{
		WinActivate
	}
	Sleep 500
	;change directory to where we want the data saved
	SendInput cd{space}'%A_LoopFileFullPath%'
	SendInput {Enter}
	;Sleep 200
	;run the save files script which at the end should exit matlab
	;must convert the save script to a function with a patient ID parameter
	;must make sure the save script is in the matlab path
	SendInput id='%A_LoopFileName%'
	SendInput {Enter}
	Sleep 200
	SendInput fSHHSExp
	SendInput {Enter}
	;need to wait until saving is done
	;Exit from Matlab, exit is built into function call
	WinWaitClose, MATLAB,,MatLab Export
	;must make sure profusion window active
	;Must click close in profusion
	IfWinExist MatLab Export
	{
		WinActivate
	}
	Click 450, 300
	Sleep 300
	SendInput {Enter}
	Sleep 1000
	;Close Profusion
	SendInput !S
	SendInput x
	WinWaitClose, Profusion
	Sleep 1000
}
return

