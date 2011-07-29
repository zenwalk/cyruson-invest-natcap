; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "InVEST"
!define PRODUCT_VERSION "2.2 beta"
!define PRODUCT_PUBLISHER "The Natural Capital Project"
!define PRODUCT_WEB_SITE "http://www.naturalcapitalproject.org"
!define MUI_HEADERIMAGE_BITMAP "natcap_logo.bmp"
SetCompressor /FINAL /SOLID lzma
SetCompressorDictSize 64


; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "natcap_logo.ico"

; Welcome page
!insertmacro MUI_PAGE_WELCOME
; License page
!insertmacro MUI_PAGE_LICENSE "license.rtf"
; Components page
!insertmacro MUI_PAGE_DIRECTORY
; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
;!define MUI_FINISHPAGE_SHOWREADME "InVEST_Documentation_v2.2.pdf"
;!insertmacro MUI_PAGE_FINISH

; Language files
!insertmacro MUI_LANGUAGE "English"

; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "InVEST_2.2_beta-Setup.exe"
InstallDir "C:\InVEST22"
ShowInstDetails show

Function .onInit
 System::Call 'kernel32::CreateMutexA(i 0, i 0, t "myMutex") i .r1 ?e'
 Pop $R0

 StrCmp $R0 0 +3
   MessageBox MB_OK|MB_ICONEXCLAMATION "The installer is already running."
   Abort

FunctionEnd


Section "InVEST Tools" SEC01
  SetOutPath "$INSTDIR\python"
  SetOverwrite try
  File /r ..\python
  
  SetOutPath "$INSTDIR"
  SetOverwrite try
  File /r ..\invest-data-tmp\*
  
SectionEnd
