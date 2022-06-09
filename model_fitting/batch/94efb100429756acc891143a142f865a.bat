@echo off
REM automatically generated
ECHO generated on host: WPIA-DIDE304
ECHO generated on date: 2022-06-09
ECHO didehpc version: 0.3.5
ECHO context version: 0.3.0
ECHO running on: %COMPUTERNAME%
set CONTEXT_WORKDRIVE=Q:
set CONTEXT_WORKDIR=COV_Italy_multistrain\model_fitting
set CONTEXT_ROOT=Q:\COV_Italy_multistrain\model_fitting
set CONTEXT_ID=48c117331d7b5e5392ba9f84cfb4127b
set R_LIBS_USER=Q:\COV_Italy_multistrain\model_fitting\lib\windows\4.1
call setr64_4_1_2.bat
REM If Java is wanted, then call setJava64.
REM If called with blank, it adds default JRE.
IF 'FALSE'=='TRUE' (
  call setJava64.bat 
)
ECHO mapping Q: -^> \\fi--san03.dide.ic.ac.uk\homes\bnc19
net use Q: \\fi--san03.dide.ic.ac.uk\homes\bnc19 /y
ECHO mapping T: -^> \\fi--didef3.dide.ic.ac.uk\tmp
net use T: \\fi--didef3.dide.ic.ac.uk\tmp /y
ECHO Using Rtools at T:\Rtools\Rtools40
set PATH=T:\Rtools\Rtools40\mingw$(WIN)\bin;T:\Rtools\Rtools40\usr\bin;%PATH%
set BINPREF=T:/Rtools/Rtools40/mingw$(WIN)/bin/
ECHO This is a parallel job: will use %CPP_NUMCPUS%
set CONTEXT_CORES=%CCP_NUMCPUS%
set REDIS_HOST=11.0.0.1
set REDIS_URL=redis://11.0.0.1:6379
%CONTEXT_WORKDRIVE%
cd \%CONTEXT_WORKDIR%
ECHO working directory: %CD%
ECHO this is a single task
set CONTEXT_TASK_ID=94efb100429756acc891143a142f865a
set CONTEXT_LOGFILE=Q:\COV_Italy_multistrain\model_fitting\logs\%CONTEXT_TASK_ID%
ECHO logfile: %CONTEXT_LOGFILE%
@REM The quoting here is necessary for paths with spaces.
ECHO on
Rscript "Q:\COV_Italy_multistrain\model_fitting\bin\task_run" "%CONTEXT_ROOT%" %CONTEXT_TASK_ID% > "%CONTEXT_LOGFILE%" 2>&1
@ECHO off
%SystemDrive%
set ErrorCode=%ERRORLEVEL%
ECHO Removing mapping Q:
net use Q: /delete /y
ECHO Removing mapping T:
net use T: /delete /y
set ERRORLEVEL=%ErrorCode%
if %ERRORLEVEL% neq 0 (
  ECHO Error running task
  EXIT /b %ERRORLEVEL%
)
@ECHO Quitting
