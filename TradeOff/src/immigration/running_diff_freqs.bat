@echo off
setlocal enabledelayedexpansion

rem Set constant values for the first, second, and fourth arguments
set arg1=10
set arg2=1
set arg4=1
set arg5=20
set arg6=25

rem Loop over a range of values for the third argument
for %%C in (10 20 40 80 160 320 640) do (
    call :retry %%C
)

endlocal
pause
exit /b

:retry
rem Initialize retry counter and max retries
set retries=0
set max_retries=3

rem Print the current combination of arguments
echo Running with arguments: %arg1% %arg2% %1 %arg4% %arg5% %arg6%

:retry_loop
rem Run the Julia script with the current combination of arguments
julia .\src\immigration\assemble_to_averages.jl %arg1% %arg2% %1 %arg4% %arg5% %arg6%

rem Check the error level
if %errorlevel% neq 0 (
    rem Increment the retry counter
    set /a retries+=1
    echo Error encountered. Attempt %retries% of %max_retries%.
    
    rem Check if max retries reached
    if %retries% lss %max_retries% (
        rem Wait for a few seconds before retrying (optional)
        timeout /t 5 /nobreak >nul
        rem Retry the command
        goto retry_loop
    ) else (
        echo Maximum retries reached. Moving to next set of arguments.
    )
)
exit /b
