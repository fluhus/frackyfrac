:: Runs frcfrc on the test data and checks the results.

echo off

set frcfrc=..\..\bin\frcfrc

for %%f in (uwtd1 uwtd2) do (
    %frcfrc% -i %%f.dense -t %%f.tree -o %%f.got
    %frcfrc% -s -i %%f.sparse -t %%f.tree -o %%f.s.got
) 2> nul

for %%f in (wtd) do (
    %frcfrc% -w -i %%f.dense -t %%f.tree -o %%f.got
    %frcfrc% -w -s -i %%f.sparse -t %%f.tree -o %%f.s.got
) 2> nul

for %%f in (uwtd1 uwtd2 wtd) do (
    fc %%f.got %%f.want
    fc %%f.s.got %%f.want
)

del *.got
