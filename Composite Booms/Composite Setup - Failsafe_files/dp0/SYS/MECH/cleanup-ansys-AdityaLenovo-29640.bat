@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 5500)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 30112)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 24452)
if /i "%LOCALHOST%"=="AdityaLenovo" (taskkill /f /pid 29640)

del /F cleanup-ansys-AdityaLenovo-29640.bat
