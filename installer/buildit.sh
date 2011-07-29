#!/bin/bash

V=$1

if [ $# -eq 0 ]
then
  V=2.2
fi

wine ~/.wine/drive_c/Program\ Files/NSIS/makensis.exe /DVERSION=$V ./invest_install_build.nsi
