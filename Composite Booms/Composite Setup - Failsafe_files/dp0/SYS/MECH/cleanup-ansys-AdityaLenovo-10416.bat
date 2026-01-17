@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 28240)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 440)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 27688)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 10416)

del /F cleanup-ansys-AdityaLenovo-10416.bat
