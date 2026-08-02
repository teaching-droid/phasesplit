@echo off
rem Build phasesplit with MSVC.
rem Visual Studio is located by walking the install folders rather than through
rem vswhere: the x86 Program Files variable has brackets in its name and those
rem upset cmd's parser inside if-blocks, which is a fight not worth having.
setlocal

if not "%VSCMD_VER%"=="" goto :have_env
where cl.exe >nul 2>&1 && goto :have_env

set "VCVARS="
for /d %%d in ("%ProgramFiles%\Microsoft Visual Studio\*") do (
    for /d %%e in ("%%d\*") do (
        if exist "%%e\VC\Auxiliary\Build\vcvars64.bat" set "VCVARS=%%e\VC\Auxiliary\Build\vcvars64.bat"
    )
)
if not defined VCVARS for /d %%d in ("%ProgramFiles(x86)%\Microsoft Visual Studio\*") do (
    for /d %%e in ("%%d\*") do (
        if exist "%%e\VC\Auxiliary\Build\vcvars64.bat" set "VCVARS=%%e\VC\Auxiliary\Build\vcvars64.bat"
    )
)
if not defined VCVARS goto :no_vs

call "%VCVARS%" >nul
if errorlevel 1 goto :no_vs

:have_env
if not exist build mkdir build
pushd build

rem /O2 speed, /fp:fast because this is audio DSP.
rem No /arch: switch on purpose - the wide SIMD paths are chosen at run time,
rem so one binary still starts on a plain SSE2 machine.
cl /nologo /W4 /WX /O2 /fp:fast /GS- /std:c11 ^
   /I..\src ^
   ..\src\main.c ..\src\wav.c ..\src\fft.c ..\src\selftest.c ..\src\split.c ..\src\cpu.c ..\src\dsp.c ..\src\thread.c ..\src\fx.c ..\src\upmix.c ^
   /Fe:phasesplit.exe /Fo:.\ ^
   /link /INCREMENTAL:NO
if errorlevel 1 goto :failed

popd
echo.
echo Built build\phasesplit.exe
exit /b 0

:failed
popd
echo.
echo BUILD FAILED
exit /b 1

:no_vs
echo Visual Studio with the C++ build tools was not found.
echo Open a "Developer Command Prompt" and run this again, or install the
echo "Desktop development with C++" workload.
exit /b 1
