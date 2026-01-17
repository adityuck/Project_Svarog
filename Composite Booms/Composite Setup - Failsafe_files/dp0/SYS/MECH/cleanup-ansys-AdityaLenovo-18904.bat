@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 29052)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 3104)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 22320)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 18904)

del /F cleanup-ansys-AdityaLenovo-18904.bat
