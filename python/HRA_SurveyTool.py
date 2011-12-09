from Tkinter import *
import tkMessageBox
import tkFileDialog
import os

def StringVarListIndex(list,item):
    """
    Returns the index of a value in a StringVar list.
    """
    for i in range(len(list)):
        if list[i].get()==item.get():
            return i
    raise ValueError, "list.index(x): x not in list"

class HRA:
    """
    Graphical user interface for habitat risk ratings.
    """
    def __init__(self):
        self.tk = Tk()
##        self.tk.resizable(width=FALSE, height=FALSE)
##        self.tk.minsize(width=640, height=480)
        #self.tk.config(background="White")

        self.weightClasses=["-","0","+"]

        #validation
        self.vcmdInt = (self.tk.register(self.OnValidateInt), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.vcmdString = (self.tk.register(self.OnValidateString), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        
        #self.tk.option_add("*font","Courier")
        
        self.tk.geometry("640x480+"+str((self.tk.winfo_screenwidth()-640)/2)+"+"+str((self.tk.winfo_screenheight()-480)/2))
        self.page=None
        self.helpRoot=None
        self.helpCanvas=None

        #set mode
        self.importedData=False
        self.dataValid=False
        self.editSurvey=False
        self.editWeights=False
        self.basicMode=True
        self.showHelp=IntVar()
        self.showHelp.set(0)

        self.promptWidth=600

        self.helpPages={}
        self.LoadHelpPages()
        self.helpKey="unknown"

        self.habitats2=None
        self.stressors2=None
        
        # create menubar
        self.menubar = Menu(self.tk)
        filemenu = Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="New", command=self.New)
##        filemenu.add_command(label="New Advanced", command=self.AdvancedMode)
##        filemenu.add_command(label="Open")
##        filemenu.add_command(label="Save")
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.tk.quit)
        self.menubar.add_cascade(label="File", menu=filemenu)

        self.indexmenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Index", menu=self.indexmenu)
        self.indexmenu.add_command(label="You must first select the assessment origin.", command=None)

        self.dataIndexBuilt=False

        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_checkbutton(label="Show Help", variable=self.showHelp, command=self.OnShowHelp)
##        self.helpmenu.add_command(label="User Guide", command=self.OnUserGuide)
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)

        # display the menu
        self.tk.config(menu=self.menubar)

        self.Home()

##        #gui flow
##        self.sequence="Survey"
##        surveyDesignSequence=[]
##        surveyNamesSequence=[self.GetDataClassNames,
##                             self.GetIntensityNames,
##                             self.GetManagementNames,
##                             self.GetAreaChangeNames,
##                             self.GetStructureChangeNames,
##                             self.GetDisturbanceFrequencyNames,
##                             self.GetOverlapTimeNames,
##                             self.GetMortalityRateNames,
##                             self.GetRecruitmentPatternNames,
##                             self.GetConnectivityNames,
##                             self.GetRecoveryTimeNames]
##        
##        dataDesignSequence=[self.GetHabitatStressors,
##                            self.GetHabitatStressorsValidate]
##
##        dataValuesSequence=[self.GetHabitatNames,
##                            self.GetStressorNames,
##                            self.GetIntensity,
##                            self.GetManagement,
##                            self.GetDistance,
##                            self.GetAreaChange,
##                            self.GetStructureChange,
##                            self.GetDisturbanceFrequency,
##                            self.GetTimeOverlap,
##                            self.GetMortalityRate,
##                            self.GetRecruitmentPattern,
##                            self.GetConnectivity,
##                            self.GetRecoveryTime]
##        
##        surveyhelpSequence=[]
##        self.sequences={"Survey":dataValuesSequence}
##        self.pageIndex=0
        
        self.tk.mainloop()

    def New(self):
        if tkMessageBox.askyesno("New", "Discard changes?", default=tkMessageBox.NO, icon=tkMessageBox.WARNING):
            self.Home()

    def LoadHelpPages(self):
        folderName = "\\".join(sys.argv[0].split("\\")[:-1])+"\\"
        fileName = "help.rst"
        helpfile = open(folderName+fileName,'r')
        rows=[row.strip().split(":\n") for row in helpfile.read().strip().split(".. _my-")]
        rows.pop(0)
        for id,text in rows:
            self.helpPages[id]=text.strip()
        helpfile.close()
    

    def OnValidateInt(self, d, i, P, s, S, v, V, W):
        """
        Validate strings to only allow integers.
        """
        try:
            int(S)
            if i=="0" and S == "0":
                return False
            else:
                return True
        except:
            return False

    def OnValidateString(self, d, i, P, s, S, v, V, W):
        """
        Validate strings to disallow commas.
        """
        return not S == ","

##    def Next(self):
##        """
##        Calls the next page in the sequence.
##        """
##        if self.pageIndex+1 < len(self.sequences[self.sequence]):
##            self.pageIndex+=1
##            self.sequences[self.sequence][self.pageIndex]()
##        else:
##            self.Home()
##
##    def NextCycle(self):
##        """
##        Calls the next page in the cycle sequence.
##        """
##        if self.cycleIndex+1 < self.cycleMax:
##            self.cycleIndex+=1
##            self.sequences[self.sequence][self.pageIndex]()
##        elif self.pageIndex+1 < len(self.sequences[self.sequence]):
##            self.pageIndex+=1
##            self.sequences[self.sequence][self.pageIndex]()
##        else:
##            self.Home()
##
##    def Previous(self):
##        """
##        Calls the previous page in the sequence
##        """
##        if self.pageIndex > 0:
##            self.pageIndex-=1
##            self.sequences[self.sequence][self.pageIndex]()
##        else:
##            self.Home()
##
##    def PreviousCycle(self):
##        """
##        Calls the previous page in the sequence
##        """
##        if self.cycleIndex > 0:
##            self.cycleIndex-=1
##            self.sequences[self.sequence][self.pageIndex]()
##        elif self.pageIndex > 0:
##            self.pageIndex-=1
##            self.sequences[self.sequence][self.pageIndex]()
##        else:
##            self.Home()
##
##    
    def OnShowHelp(self):
        """
        Toggeles the show help state.
        """
        menuBarHeight=45
        helpWidth=220
        self.tk.update_idletasks()
        
        if self.showHelp.get():
##            if not self.helpRoot==None:
##                self.helpRoot.deiconify()
##                self.helpRoot.geometry("220x480+"+str(self.tk.winfo_rootx()+self.tk.winfo_width())+"+"+str(self.tk.winfo_rooty()-menuBarHeight))  
##            else:
            self.helpRoot = Toplevel(self.tk)
            self.helpRoot.wm_title("Help")

            #position help window next to left of wizard matching height
            width,height,xoffset,yoffset=map(int,self.tk.winfo_geometry().replace("x","+").split("+"))
            
            self.helpRoot.geometry(str(helpWidth)+"x"+str(height)+"+"+str(xoffset+width+5)+"+"+str(yoffset))
            menubar = Menu(self.helpRoot)
            self.helpRoot.config(menu=menubar)
            self.helpRoot.resizable(width=FALSE, height=TRUE)
            self.helpRoot.protocol("WM_DELETE_WINDOW",self.OnHideHelp)

            #create canvas with scroll bar
            vscrollbar = Scrollbar(self.helpRoot, orient=VERTICAL)
            vscrollbar.grid(row=0, column=1, sticky=N+S)
            self.helpCanvas = Canvas(self.helpRoot, yscrollcommand=vscrollbar.set)
            self.helpCanvas.grid(row=0, column=0, sticky=N+S+E+W)
            #self.helpCanvas.config(background="White")
            vscrollbar.config(command=self.helpCanvas.yview)

            #expand canvas
            self.helpRoot.grid_rowconfigure(0, weight=1)
            self.helpRoot.grid_columnconfigure(0, weight=1)

            # create canvas contents

            frame = Frame(self.helpCanvas)
            frame.rowconfigure(1, weight=1)
            frame.columnconfigure(1, weight=1)

##                rows = 30
##                for i in range(2,rows):
##                    for j in range(1,10):
##                        button = Button(frame, text="[%d,%d]" % (i,j))
##                        button.grid(row=i, column=j, sticky='news')

            #frame.pack.forget()
            Label(frame, text=self.helpPages[self.helpKey], wraplength=helpWidth-20, justify=LEFT).grid(row=0, column=0, sticky=W)

            self.helpCanvas.create_window(0, 0, anchor=NW, window=frame)

            frame.update_idletasks()

            self.helpCanvas.config(scrollregion=self.helpCanvas.bbox("all"))

##            button = Button(frame, text="Hi")
##            button.grid(row=0, column=0)
        elif not self.helpRoot == None:
            self.helpRoot.withdraw()


    def OnHideHelp(self):
        """
        Unchecks the show help when the help window is closed.
        """
        self.showHelp.set(False)
        self.OnShowHelp()

    def OnUserGuide(self):
        os.system("start C:\\InVEST_2.1.0\\InVEST_2.1.0_Documentation.pdf")

##    def BasicMode(self):
##        self.basicMode=True
##        self.Home()
##
##    def AdvancedMode(self):
##        self.basicMode=False
##        self.Home()

    def Home(self):
        """
        Calls the basic or advanced menu based on the user mode.
        """
        self.tk.wm_title("Marine InVEST - Habitat Risk Assesment")
        self.menubar.delete(2)
        self.dataIndexBuilt=False

        self.indexmenu = Menu(self.menubar, tearoff=0)
        self.menubar.insert_cascade(2, label="Index", menu=self.indexmenu)
        self.indexmenu.add_command(label="You must first select the assessment origin.", command=None)

        if self.basicMode:
            self.BasicMenu()
        else:
            self.AdvancedMenu()

            
    def NewData(self):
        """
        Initializes new variables for data.
        """
        self.importedData=False
        self.habitats=StringVar()
        self.stressors=StringVar()
        self.DefaultWeights()

        self.indexmenu.delete(1)        
        self.criteriaIndex()
        self.indexmenu.add_command(label="The complete index will appear after entering habitat and stressor names.", command=None)

        self.DefaultEditSurvey()
   
    def NewSurvey(self):
        """
        Initializes new variables for a survey.
        """
        self.dataClassNames=[]
        self.dataClasses=0
        
        self.stressorIntensityNames=[]
        self.stressorIntensities=0
        
        self.stressorManagementNames=[]
        self.stressorManagements=0
        
        self.areaChangeNames=[]
        self.areaChanges=0
        
        self.structureChangeNames=[]
        self.structureChanges=0
        
        self.disturbanceFrequencyNames=[]
        self.disturbanceFrequencies=0
        
        self.overlapTimeNames=[]
        self.overlapTimes=0
        
        self.mortalityRateNames=[]
        self.mortalityRates=0
        
        self.recruitmentPatternNames=[]
        self.recruitmentPatterns=0
        
        self.connectivityNames=[]
        self.connectivities=0
        
        self.recoveryTimeNames=[]
        self.recoveryTimes=0


    #example data
    def DefaultSurvey(self):
        """
        Sets the survey to the WCVI example.
        """
        dataClassNames=["unknown","best", "adequate", "limited"]      
        stressorIntensityNames=["no score","low","medium","high"]
        stressorManagementNames=["no score","very effective", "somewhat effective", "not effective"]
        areaChangeNames=["no score","low (0-20%)","medium (20-50%)","high (50-100%)"]
        structureChangeNames=["no score","low (0-20%)","medium (20-50%)","high (50-100%)"]
        disturbanceFrequencyNames=["no score","daily to weekly","several times per year","annually or less often"]
        overlapTimeNames=["no score","0-4 months","4-8 months","8-12 months"]
        mortalityRateNames=["no score", "high (>=80%)", "moderate (20-50%)", "low (0-20%)"]
        recruitmentPatternNames=["no score","annual or more often","every 1-2 years","every 2+ years"]
        connectivityNames=["no score","high dispersal (100km+)", "medium dispersal (10-100km)", "low dispersal (<10km)"]
        recoveryTimeNames=["no score","0-1 years","1-10 years","10+ years"]
        
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,
                        stressorManagementNames,
                        areaChangeNames,
                        structureChangeNames,
                        disturbanceFrequencyNames,
                        overlapTimeNames,
                        mortalityRateNames,
                        recruitmentPatternNames,
                        connectivityNames,
                        recoveryTimeNames)

    def DefaultData(self):
        """
        Sets the data to the WCVI example.
        """
        self.DefaultSurvey()

        assessmentName="WCVI"
        
        habitatNames=["kelp","eelgrass","hard bottom","soft bottom"]
        habitatQuality=[1,1,2,2]
                  
        stressorNames=["Finfish Aquaculture","Shellfish Aquaculture","Comm Salmon (Troll)","Rec Fishing"]
        stressorQuality=[1,1,2,2]
                  
        stressorIntensity=[3,2,3,2]
        stressorIntensityQuality=[1,1,1,1]

        stressorManagement=[3,2,3,3]
        stressorManagementQuality=[1,2,2,2]

        stressorDistance=[300,250,150,100]

        areaChange=[2,2,1,1,3,3,2,2,0,0,0,0,0,0,0,0]
        areaChangeQuality=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]

        structureChange=[2,2,1,1,3,3,2,2,1,1,1,1,1,1,1,1]
        structureChangeQuality=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]

        disturbanceFrequency=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        disturbanceFrequencyQuality=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]

        overlapTime=[3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2]
        overlapTimeQuality=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

        mortalityRate=[1,1,0,0]
        mortalityRateQuality=[1,1,0,0]

        recruitmentPattern=[1,1,0,0]
        recruitmentPatternQuality=[1,1,1,1]

        connectivity=[2,1,0,0]
        connectivityQuality=[1,1,0,0]

        recoveryTime=[1,1,3,1]
        recoveryTimeQuality=[1,1,1,1]


        self.LoadData(assessmentName,
                      habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)


    def DefaultWeights(self):
        weights = [1,1,1,1,1,1,1,1,1,1,1]
        self.LoadWeights(weights)

    def LoadWeights(self,weights):
        self.weights=self.IntVarList(weights)
    
    def WCVI(self):
        """
        Calls DefaultSurvey and DefaultData.
        """
        self.DefaultSurvey()
        self.DefaultData()
        self.DefaultWeights()

        #build index
        self.menubar.delete(2)        
        self.menubar.insert_cascade(2,label="Index",menu=self.indexmenu)
        self.indexmenu.delete(0)
        self.criteriaIndex()
        self.dataIndex()
        self.weightIndex()

    def criteriaIndex(self):
        criteria = Menu(self.indexmenu, tearoff=0)
        criteria.add_command(label="Data Quality", command=self.GetDataClassNames)
        criteria.add_command(label="Intensity", command=self.GetIntensityNames)
        criteria.add_command(label="Management", command=self.GetManagementNames)
        criteria.add_command(label="Area", command=self.GetAreaChangeNames)
        criteria.add_command(label="Structure", command=self.GetStructureChangeNames)
        criteria.add_command(label="Disturbance", command=self.GetDisturbanceFrequencyNames)
        criteria.add_command(label="Temporal", command=self.GetOverlapTimeNames)
        criteria.add_command(label="Mortality", command=self.GetMortalityRateNames)
        criteria.add_command(label="Recruitment", command=self.GetRecruitmentPatternNames)
        criteria.add_command(label="Connectivity", command=self.GetConnectivityNames)
        criteria.add_command(label="Recovery", command=self.GetRecoveryTimeNames)                     
        self.indexmenu.add_cascade(label="Criteria", menu=criteria)
        
    def dataIndex(self):
        self.dataIndexBuilt=True
        data = Menu(self.indexmenu, tearoff=0)
        data.add_command(label="Quantity", command=self.GetHabitatStressors)
        data.add_command(label="Habitats", command=self.GetHabitatNames)
        data.add_command(label="Stressors", command=self.GetStressorNames)
        data.add_command(label="Intensity", command=self.GetIntensity)
        data.add_command(label="Management", command=self.GetManagement)
        data.add_command(label="Distance", command=self.GetDistance)

        habitatRange = range(int(self.habitats.get()))
        for label,function in [("Area Change",self.GetHabitatAreaChangeIndex),
                               ("Structure Change",self.GetHabitatStructureChangeIndex),
                               ("Disturbance Frequencey",self.GetHabitatDisturbanceFrequencyIndex),
                               ("Temporal Overlap",self.GetHabitatTimeOverlapIndex)]:
            tmp = Menu(data, tearoff=0)
            for i in habitatRange:
                tmp.add_command(label=self.habitatNames[i].get(), command=lambda i=i: function(i))
            data.add_cascade(label=label, menu=tmp)

        data.add_command(label="Mortality", command=self.GetMortalityRate)
        data.add_command(label="Recruitment", command=self.GetRecruitmentPattern)
        data.add_command(label="Connectivity", command=self.GetConnectivity)
        data.add_command(label="Maturity/Recovery", command=self.GetRecoveryTime)
        
        self.indexmenu.add_cascade(label="Data", menu=data)

    def weightIndex(self):
        weighting = Menu(self.indexmenu, tearoff=0)
        weighting.add_command(label="Exposure", command=self.GetExposureWeights)
        weighting.add_command(label="Consequence", command=self.GetConsequenceWeights)
        self.indexmenu.add_cascade(label="Weighting", menu=weighting)

    #I/O
    def ReadFile(self,filename):
        """
        Reads data in the CSV format.
        """
        infile=open(filename,"r")
        assessmentName,survey,habitat,stressor,habitatStressor,weights=infile.read().strip().split("\n\n")
        infile.close()

        #survey
        survey=survey.split("\n")
        row=survey.pop(0).split(",")
        dataClasses=int(row.pop(0))
        dataClassNames=row

        row=survey.pop(0).split(",")
        stressorIntensities=int(row.pop(0))
        stressorIntensityNames=row

        row=survey.pop(0).split(",")
        stressorManagements=int(row.pop(0))
        stressorManagementNames=row

        row=survey.pop(0).split(",")
        areaChanges=int(row.pop(0))
        areaChangeNames=row

        row=survey.pop(0).split(",")
        structureChanges=int(row.pop(0))
        structureChangeNames=row

        row=survey.pop(0).split(",")
        disturbanceFrequencies=int(row.pop(0))
        disturbanceFrequencyNames=row

        row=survey.pop(0).split(",")
        overlapTimes=int(row.pop(0))
        overlapTimeNames=row
        
        row=survey.pop(0).split(",")
        mortalityRates=int(row.pop(0))
        mortalityRateNames=row

        row=survey.pop(0).split(",")
        recruitmentPatterns=int(row.pop(0))
        recruitmentPatternNames=row
        
        row=survey.pop(0).split(",")
        connectivities=int(row.pop(0))
        connectivityNames=row

        row=survey.pop(0).split(",")
        recoveryTimes=int(row.pop(0))
        recoveryTimeNames=row
        
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,
                        stressorManagementNames,
                        areaChangeNames,
                        structureChangeNames,
                        disturbanceFrequencyNames,
                        overlapTimeNames,
                        mortalityRateNames,
                        recruitmentPatternNames,
                        connectivityNames,
                        recoveryTimeNames)


        #habitat data
        habitat=habitat.split("\n")
        habitat.pop(0) #skip number of habitats
        habitat.pop(0) #skip column headers

        habitatNames=[]
        habitatQuality=[]
        mortalityRate=[]
        mortalityRateQuality=[]
        recruitmentPattern=[]
        recruitmentPatternQuality=[]
        connectivity=[]
        connectivityQuality=[]
        recoveryTime=[]
        recoveryTimeQuality=[]
        for row in habitat:
            row=row.split(",")
            habitatNames.append(row[1])
            habitatQuality.append(int(row[2]))
            mortalityRate.append(int(row[3]))
            mortalityRateQuality.append(int(row[4]))
            recruitmentPattern.append(int(row[5]))
            recruitmentPatternQuality.append(int(row[6]))
            connectivity.append(int(row[7]))
            connectivityQuality.append(int(row[8]))
            recoveryTime.append(int(row[9]))
            recoveryTimeQuality.append(int(row[10]))
        
        #stressor data
        stressor=stressor.split("\n")
        stressor.pop(0) #skip number of stressors
        stressor.pop(0) #skip header row

        stressorNames=[]
        stressorQuality=[]       
        stressorIntensity=[]
        stressorIntensityQuality=[]
        stressorManagement=[]
        stressorManagementQuality=[]
        stressorDistance=[]
        for row in stressor:
            row=row.split(",")
            stressorNames.append(row[1])
            stressorQuality.append(int(row[2]))
            stressorIntensity.append(int(row[3]))
            stressorIntensityQuality.append(int(row[4]))
            stressorManagement.append(int(row[5]))
            stressorManagementQuality.append(int(row[6]))
            stressorDistance.append(float(row[7]))

        #habitat-stressor data
        habitatStressor=habitatStressor.split("\n")
        habitatStressor.pop(0) #skip number or habitats*stressors
        habitatStressor.pop(0) #skip header

        areaChange=[]
        areaChangeQuality=[]
        structureChange=[]
        structureChangeQuality=[]
        disturbanceFrequency=[]
        disturbanceFrequencyQuality=[]
        overlapTime=[]
        overlapTimeQuality=[]
        for row in habitatStressor:
            row=row.split(",")
            areaChange.append(int(row[4]))
            areaChangeQuality.append(int(row[5]))
            structureChange.append(int(row[6]))
            structureChangeQuality.append(int(row[7]))
            disturbanceFrequency.append(int(row[8]))
            disturbanceFrequencyQuality.append(int(row[9]))
            overlapTime.append(int(row[10]))
            overlapTimeQuality.append(int(row[11]))

        self.LoadData(assessmentName,
                      habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)

        #weights
        self.LoadWeights(map(int,weights.split("\n")))
        

    def GetVarValue(self,variable):
        """
        Calls the .get() function on the parameter.
        """
        return variable.get()
            
    def WriteFile(self,filename):
        """
        Writes data to the CSV format.
        """
        outfile=open(filename,"w")

        #assessment name
        outfile.write(self.assessmentName.get()+"\n\n")

        #survey        
        outfile.write(self.dataClasses.get()+","+",".join(map(self.GetVarValue, self.dataClassNames)))
        outfile.write("\n"+self.stressorIntensities.get()+","+",".join(map(self.GetVarValue, self.stressorIntensityNames)))
        outfile.write("\n"+self.stressorManagements.get()+","+",".join(map(self.GetVarValue, self.stressorManagementNames)))
        outfile.write("\n"+self.areaChanges.get()+","+",".join(map(self.GetVarValue, self.areaChangeNames)))
        outfile.write("\n"+self.structureChanges.get()+","+",".join(map(self.GetVarValue, self.structureChangeNames)))
        outfile.write("\n"+self.disturbanceFrequencies.get()+","+",".join(map(self.GetVarValue, self.disturbanceFrequencyNames)))
        outfile.write("\n"+self.overlapTimes.get()+","+",".join(map(self.GetVarValue, self.overlapTimeNames)))
        outfile.write("\n"+self.mortalityRates.get()+","+",".join(map(self.GetVarValue, self.mortalityRateNames)))
        outfile.write("\n"+self.recruitmentPatterns.get()+","+",".join(map(self.GetVarValue, self.recruitmentPatternNames)))
        outfile.write("\n"+self.connectivities.get()+","+",".join(map(self.GetVarValue, self.connectivityNames)))
        outfile.write("\n"+self.recoveryTimes.get()+","+",".join(map(self.GetVarValue, self.recoveryTimeNames)))
        outfile.write("\n\n")

         #habitat data
        outfile.write(self.habitats.get()+"\n")
        outfile.write("Habitat ID,Habitat Name,Habitat DQ,Mortality,Mortality DQ,Recruitment,Recruitment DQ,Connectivity,Connectivity DQ,Regeneration,Regeneration DQ")
        for i in range(int(self.habitats.get())):
            outfile.write("\n"+str(i+0))
            outfile.write(","+self.habitatNames[i].get())
            outfile.write(","+str(self.habitatQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.mortalityRateNames,self.mortalityRate[i])))
            outfile.write(","+str(self.mortalityRateQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.recruitmentPatternNames,self.recruitmentPattern[i])))
            outfile.write(","+str(self.recruitmentPatternQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.connectivityNames,self.connectivity[i])))
            outfile.write(","+str(self.connectivityQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.recoveryTimeNames,self.recoveryTime[i])))
            outfile.write(","+str(self.recoveryTimeQuality[i].get()+0))
        outfile.write("\n\n")
        
        #stressor data
        outfile.write(self.stressors.get()+"\n")
        outfile.write("Stressor ID,Stressor Name,Stressor DQ,Intensity,Intensity DQ,Management,Managment DQ,Buffer")
        for i in range(int(self.stressors.get())):
            outfile.write("\n"+str(i+0))
            outfile.write(","+self.stressorNames[i].get())
            outfile.write(","+str(self.stressorQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.stressorIntensityNames,self.stressorIntensity[i])))
            outfile.write(","+str(self.stressorIntensityQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.stressorManagementNames,self.stressorManagement[i])))
            outfile.write(","+str(self.stressorManagementQuality[i].get()+0))
            outfile.write(","+str(self.stressorDistance[i].get()))
        outfile.write("\n\n")

        #habitat-stressor data
        outfile.write(str(int(self.habitats.get())*int(self.stressors.get()))+"\n")
        outfile.write("Habitat ID,Habitat Name,Stressor ID,Stressor Name,Area Change,Area Change DQ,Structure Change,Structure Change DQ,Disturbance Frequency,Disturbance Frequency DQ,Temporal Overlap,Temporal Overlap DQ")
        for i in range(int(self.habitats.get())):
            for j in range(int(self.stressors.get())):
                outfile.write("\n"+str(i+0))
                outfile.write(","+self.habitatNames[i].get())
                outfile.write(","+str(j+0))
                outfile.write(","+self.stressorNames[j].get())
                outfile.write(","+str(StringVarListIndex(self.areaChangeNames,self.areaChange[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.areaChangeQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.structureChangeNames,self.structureChange[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.structureChangeQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.disturbanceFrequencyNames,self.disturbanceFrequency[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.disturbanceFrequencyQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.overlapTimeNames,self.overlapTime[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.overlapTimeQuality[i*int(self.stressors.get())+j].get()+0))
        outfile.write("\n\n")

        #weights data
        outfile.write("\n".join(map(str,map(self.GetVarValue, self.weights))))
        
        outfile.close()

    def LoadValues(self,dataClassNames,
             stressorIntensityNames,
             stressorManagementNames,
             areaChangeNames,
             structureChangeNames,
             disturbanceFrequencyNames,
             overlapTimeNames,
             mortalityRateNames,
             recruitmentPatternNames,
             connectivityNames,
             recoveryTimeNames,
             assessmentName,
             habitatNames,habitatQuality,
             stressorNames,stressorQuality,
             stressorIntensity,stressorIntensityQuality,
             stressorManagement,stressorManagementQuality,
             stressorDistance,
             areaChange,areaChangeQuality,
             structureChange,structureChangeQuality,
             disturbanceFrequency,disturbanceFrequencyQuality,
             overlapTime,overlapTimeQuality,
             mortalityRate,mortalityRateQuality,
             recruitmentPattern,recruitmentPatternQuality,
             connectivity,connectivityQuality,
             recoveryTime,recoveryTimeQuality):
        """
        Sets the survey and data to the specified values.
        """
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,stressorManagementNames,
                        areaChangeNames,structureChangeNames,
                        disturbanceFrequencyNames,overlapTimeNames,
                        mortalityRateNames,recruitmentPatternNames,
                        connectivityNames,recoveryTimeNames)
        self.LoadData(assessmentName,
                      habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)

    def LoadSurvey(self,dataClassNames,
                   stressorIntensityNames,
                   stressorManagementNames,
                   areaChangeNames,
                   structureChangeNames,
                   disturbanceFrequencyNames,
                   overlapTimeNames,
                   mortalityRateNames,
                   recruitmentPatternNames,
                   connectivityNames,
                   recoveryTimeNames):
        """
        Creates the data structures for a survey and sets the specified values.
        """
        #survey        
        self.dataClasses=self.StringVar(len(dataClassNames))
        self.dataClassNames=self.StringVarList(dataClassNames)
        self.stressorIntensities=self.StringVar(len(stressorIntensityNames))
        self.stressorIntensityNames=self.StringVarList(stressorIntensityNames)
        self.stressorManagements=self.StringVar(len(stressorManagementNames))
        self.stressorManagementNames=self.StringVarList(stressorManagementNames)
        self.areaChanges=self.StringVar(len(areaChangeNames))
        self.areaChangeNames=self.StringVarList(areaChangeNames)
        self.structureChanges=self.StringVar(len(structureChangeNames))
        self.structureChangeNames=self.StringVarList(structureChangeNames)
        self.disturbanceFrequencies=self.StringVar(len(disturbanceFrequencyNames))
        self.disturbanceFrequencyNames=self.StringVarList(disturbanceFrequencyNames)
        self.overlapTimes=self.StringVar(len(overlapTimeNames))
        self.overlapTimeNames=self.StringVarList(overlapTimeNames)
        self.mortalityRates=self.StringVar(len(mortalityRateNames))
        self.mortalityRateNames=self.StringVarList(mortalityRateNames)
        self.recruitmentPatterns=self.StringVar(len(recruitmentPatternNames))
        self.recruitmentPatternNames=self.StringVarList(recruitmentPatternNames)
        self.connectivities=self.StringVar(len(connectivityNames))
        self.connectivityNames=self.StringVarList(connectivityNames)
        self.recoveryTimes=self.StringVar(len(recoveryTimeNames))
        self.recoveryTimeNames=self.StringVarList(recoveryTimeNames)

    def StringVar(self,value):
        """
        Returns a StringVar variable with the specified value.
        """
        temp=StringVar()
        temp.set(str(value))
        return temp

    def StringVarList(self,valueList):
        """
        Returns a list of StringVar variables with the specified values.
        """
        return self.VarList(StringVar,valueList)

    def IntVarList(self,valueList):
        """
        Returns a list of IntVar variables with the specified values.
        """
        return self.VarList(IntVar,valueList)

    def DoubleVarList(self,valueList):
        """
        Returns a list of DoubleVar variables with the specified values.
        """
        return self.VarList(DoubleVar,valueList)

    def VarList(self,function,valueList):
        """
        Returns a list of the specfied variable type with the specified values.
        """
        varList=[]
        for i,value in enumerate(valueList):
            varList.append(function())
            varList[i].set(value)
        return varList

    def NameVarList(self,valueList,nameList):
        """
        Creates a list of StringVar variables with names specified by index.
        """
        varList=[]
        for i,value in enumerate(valueList):
            varList.append(StringVar())
            varList[i].set(nameList[value].get())
        return varList

    def LoadData(self,assessmentName,
                 habitatNames,habitatQuality,
                 stressorNames,stressorQuality,
                 stressorIntensity,stressorIntensityQuality,
                 stressorManagement,stressorManagementQuality,
                 stressorDistance,
                 areaChange,areaChangeQuality,
                 structureChange,structureChangeQuality,
                 disturbanceFrequency,disturbanceFrequencyQuality,
                 overlapTime,overlapTimeQuality,
                 mortalityRate,mortalityRateQuality,
                 recruitmentPattern,recruitmentPatternQuality,
                 connectivity,connectivityQuality,
                 recoveryTime,recoveryTimeQuality):
        """
        Creates the data structures for survey data and sets the specified values.
        """
        self.assessmentName = self.StringVar(assessmentName)
        self.habitats = self.StringVar(len(habitatNames))
        self.stressors = self.StringVar(len(stressorNames))
        self.habitatNames = self.StringVarList(habitatNames)
        self.habitatQuality = self.IntVarList(habitatQuality)
        self.stressorNames = self.StringVarList(stressorNames)
        self.stressorQuality = self.IntVarList(stressorQuality)
        self.stressorIntensity = self.NameVarList(stressorIntensity,self.stressorIntensityNames)
        self.stressorIntensityQuality = self.IntVarList(stressorIntensityQuality)
        self.stressorManagement = self.NameVarList(stressorManagement,self.stressorManagementNames)
        self.stressorManagementQuality = self.IntVarList(stressorManagementQuality)
        self.stressorDistance = self.DoubleVarList(stressorDistance)
        self.areaChange = self.NameVarList(areaChange,self.areaChangeNames)    
        self.areaChangeQuality = self.IntVarList(areaChangeQuality)
        self.structureChange = self.NameVarList(structureChange,self.structureChangeNames)
        self.structureChangeQuality = self.IntVarList(structureChangeQuality)             
        self.disturbanceFrequency = self.NameVarList(disturbanceFrequency,self.disturbanceFrequencyNames)
        self.disturbanceFrequencyQuality = self.IntVarList(disturbanceFrequencyQuality)
        self.overlapTime = self.NameVarList(overlapTime,self.overlapTimeNames)
        self.overlapTimeQuality = self.IntVarList(overlapTimeQuality)
        self.mortalityRate = self.NameVarList(mortalityRate,self.mortalityRateNames)
        self.mortalityRateQuality = self.IntVarList(mortalityRateQuality)
        self.recruitmentPattern = self.NameVarList(recruitmentPattern,self.recruitmentPatternNames)
        self.recruitmentPatternQuality = self.IntVarList(recruitmentPatternQuality)
        self.connectivity = self.NameVarList(connectivity,self.connectivityNames)
        self.connectivityQuality = self.IntVarList(connectivityQuality)
        self.recoveryTime = self.NameVarList(recoveryTime,self.recoveryTimeNames)
        self.recoveryTimeQuality = self.IntVarList(recoveryTimeQuality)
  
    
##    def LoadWeights(self,intensityWeight,managementWeight,areaWeight,structureWeight,disturbanceWeight,timeWeight,mortalityWeight,recruitmentWeight,connectivityWeight,recoveryWeight):
##        self.intensityWeight=intensityWeight
##        self.managementWeight=managementWeight
##        self.areaWeight=areaWeight
##        self.structureWeight=structureWeight
##        self.disturbanceWeight=disturbanceWeight
##        self.timeWeight=timeWeight
##        self.mortalityWeight=mortalityWeight
##        self.recruitmentWeight=recruitmentWeight
##        self.connectivityWeight=connectivityWeight
##        self.recoveryWeight=recoveryWeight
        
    #GUI
    def BasicMenu(self):
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        Label(self.page, text="How do you want to score criteria for the risk assessment?").grid(row=0,column=0,columnspan=5)

        Button(self.page, text="Create new scores", command=self.GetName, width=22).grid(row=1,column=1,sticky=S)
        Button(self.page, text="Import existing scores", command=self.Import, width=22).grid(row=1,column=2,sticky=S)
        Button(self.page, text="Sample scores", command=self.DefaultAssessment, width=22).grid(row=1,column=3,sticky=S)

        self.page.pack()

        self.helpKey="menu-basic"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetName(self):
        self.importedData=False
        
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        self.assessmentName=StringVar()
        Label(self.page, text="Enter a file name for the output (e.g. Scores_WCVI_scenarioA)").grid(row=0,column=0,columnspan=3)
        Entry(self.page, width=20, textvariable=self.assessmentName, validate="key", validatecommand=self.vcmdString).grid(row=1,column=1)

        #controls
        Button(self.page, text="Back", command=self.Home, width=4).grid(row=2,column=0,sticky=SW)
        Button(self.page, text="Next", command=self.NewData, width=4).grid(row=2,column=2,sticky=SE)

        self.page.pack()

        self.helpKey="assessment-name"
        if self.showHelp.get():
            self.OnHideHelp()
            self.showHelp.set(True)
            self.OnShowHelp()
            
##    def BasicMenu(self):
##        if not self.page == None:
##            self.page.pack_forget()
##
##        self.page = Frame(self.tk)
##
##        Label(self.page, text="Choose your survey:").grid(row=0,column=0,columnspan=4)
##        Button(self.page, text="Default", command=self.DefaultEditSurvey, width=7).grid(row=1,column=1,sticky=S)
##        Button(self.page, text="Import", command=self.Import, width=7).grid(row=1,column=2,sticky=S)
##
##        self.page.pack()
##
##        if self.showHelp.get():
##            frame = Frame(self.helpCanvas)
##            frame.rowconfigure(1, weight=1)
##            frame.columnconfigure(1, weight=1)
##
##            rows = 30
##            for i in range(2,rows):
##                for j in range(1,10):
##                    button = Button(frame, text="[%d,%d]" % (i,j))
##                    button.grid(row=i, column=j, sticky='news')
##
##            self.helpCanvas.create_window(0, 0, anchor=NW, window=frame)
##
##            frame.update_idletasks()
##
##            self.helpCanvas.config(scrollregion=self.helpCanvas.bbox("all"))

    def DefaultAssessment(self):
        self.importedData=True
        self.WCVI()
        self.EditSurvey()

    def DefaultEditSurvey(self):
        """
        Loads default survey and calls EditSurvey().
        """
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get())
        self.DefaultSurvey()
        self.EditSurvey()
    
##    def AdvancedMenu(self):
##        if not self.page == None:
##            self.page.pack_forget()
##
##        self.page = Frame(self.tk)
##
##        Label(self.page, text="Survey", width=7).grid(row=0,column=0)        
##        Button(self.page, text="New", command=self.GetNewSurvey, width=7).grid(row=1,column=0)
##        Button(self.page, text="Save", width=7).grid(row=2,column=0)
##        Button(self.page, text="Open", width=7).grid(row=3,column=0)
##        Button(self.page, text="Default", command=self.DefaultSurvey, width=7).grid(row=4,column=0)
##
##        Label(self.page, text="Data", width=7).grid(row=0,column=1)
##        Button(self.page, text="New", command=self.NewData, width=7).grid(row=1,column=1)
##        Button(self.page, text="Import", command=self.Import, width=7).grid(row=2,column=1)
##        Button(self.page, text="Export", width=7, command=self.Export).grid(row=3,column=1)
##        Button(self.page, text="Example", command=self.DefaultData, width=7).grid(row=4,column=1)
##
##        self.page.pack()

    def EditSurvey(self):
        """
        Asks user if they want to edit the survey labels.
        """
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")

        Label(self.page, text="Would you like to edit categories for scoring criteria?").grid(row=0,column=0,columnspan=4)
        Button(self.page, text="Yes", command=self.GetDataClassNames, width=3).grid(row=1,column=1,sticky=W)
        Button(self.page, text="No", command=self.GetNewData, width=3).grid(row=1,column=2,sticky=E)

        self.page.pack()            

        self.helpKey="criteria"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    #survey GUI
##    def GetNewSurvey(self):
##        self.editSurvey = False
##        self.dataClasses = StringVar()
##        self.stressorIntensities = StringVar()
##        self.stressorManagements = StringVar()
##        self.areaChanges = StringVar()
##        self.structureChanges = StringVar()
##        self.disturbanceFrequencies = StringVar()
##        self.overlapTimes = StringVar()
##        self.mortalityRates = StringVar()
##        self.recruitmentPatterns = StringVar()
##        self.connectivities = StringVar()
##        self.recoveryTimes = StringVar()
##        self.GetSurvey()

##    def GetSurvey(self):
##        self.page.pack_forget()
##            
##        self.page = Frame(self.tk)
##        Label(self.page, text="Specify the number of possible values for each survey paramter.").grid(row=0, column=0, columnspan=6)
##
##        Label(self.page, text="Data Quality").grid(row=1,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.dataClasses).grid(row=1,column=2)
##
##        Label(self.page, text="Intensity").grid(row=2,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.stressorIntensities).grid(row=2,column=2)
##
##        Label(self.page, text="Management").grid(row=2,column=3,sticky=E)
##        Entry(self.page, width=2, textvariable=self.stressorManagements).grid(row=2,column=4)
##
##        Label(self.page, text="Area Change").grid(row=3,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.areaChanges).grid(row=3,column=2)
##
##        Label(self.page, text="Structure Change").grid(row=3,column=3,sticky=E)
##        Entry(self.page, width=2, textvariable=self.structureChanges).grid(row=3,column=4)
##
##        Label(self.page, text="Disturbance Frequency").grid(row=4,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.disturbanceFrequencies).grid(row=4,column=2)
##
##        Label(self.page, text="Time Overlap").grid(row=4,column=3,sticky=E)
##        Entry(self.page, width=2, textvariable=self.overlapTimes).grid(row=4,column=4)
##
##        Label(self.page, text="Mortality Rate").grid(row=5,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.mortalityRates).grid(row=5,column=2)
##
##        Label(self.page, text="Recruitment").grid(row=5,column=3,sticky=E)
##        Entry(self.page, width=2, textvariable=self.recruitmentPatterns).grid(row=5,column=4)
##
##        Label(self.page, text="Connectivity").grid(row=6,column=1,sticky=E)
##        Entry(self.page, width=2, textvariable=self.connectivities).grid(row=6,column=2)
##
##        Label(self.page, text="Recovery Time").grid(row=6,column=3,sticky=E)
##        Entry(self.page, width=2, textvariable=self.recoveryTimes).grid(row=6,column=4)
##
##        #controls
##        Button(self.page, text="Back", command=self.Home, width=4).grid(row=7,column=0,sticky=W)
##        Button(self.page, text="Next",command=self.GetSurveyValidate, width=4).grid(row=7,column=5,sticky=E)
##
##        self.page.pack()
##
##    def StringVarN(self,n):
##        """
##        Creates a n-length list of StringVar variables.
##        """
##        temp=[]
##        for i in range(n):
##            temp.append(StringVar())
##        return temp
##
##    def GetSurveyValidate(self):
##        """
##        Validates the number of options for each survey parameter and creates the data structure to store the values.
##        """
##        try:
##            if self.editSurvey:
##                pass
##            else:
##                self.dataClassNames = self.StringVarN(int(self.dataClasses.get()))
##                self.stressorIntensityNames = self.StringVarN(int(self.stressorIntensities.get()))
##                self.stressorManagementNames = self.StringVarN(int(self.stressorManagements.get()))
##                self.areaChangeNames = self.StringVarN(int(self.areaChanges.get()))
##                self.structureChangeNames = self.StringVarN(int(self.structureChanges.get()))
##                self.disturbanceFrequencyNames = self.StringVarN(int(self.disturbanceFrequencies.get()))
##                self.overlapTimeNames = self.StringVarN(int(self.overlapTimes.get()))
##                self.mortalityRateNames = self.StringVarN(int(self.mortalityRates.get()))
##                self.recruitmentPatternNames = self.StringVarN(int(self.recruitmentPatterns.get()))
##                self.connectivityNames = self.StringVarN(int(self.connectivities.get()))
##                self.recoveryTimeNames = self.StringVarN(int(self.recoveryTimes.get()))
##                    
##                self.GetDataClassNames()
##        except ValueError:
##            tkMessageBox.showerror("Value Error","All values must be integers.")
##            self.GetSurvey()

    def GetDataClassNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'data quality\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.dataClasses.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.dataClassNames[i]).grid(row=i+1, column=2, sticky=W)

        #controls
        Button(self.page, text="Back", command=self.EditSurvey, width=4).grid(row=int(self.dataClasses.get())+1,column=0,sticky=SW)
        Button(self.page, text="Next", command=self.GetIntensityNames, width=4).grid(row=int(self.dataClasses.get())+1,column=3,sticky=SE)

        self.page.pack()
       
        self.helpKey="criteria-data"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetIntensityNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'stressor intensity\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.stressorIntensities.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.stressorIntensityNames[i]).grid(row=i+1, column=2, sticky=W)
            
        #controls
        Button(self.page, text="Back", command=self.GetDataClassNames, width=4).grid(row=int(self.stressorIntensities.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetManagementNames, width=4).grid(row=int(self.stressorIntensities.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-intensity"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
        

    def GetManagementNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'management effectiveness\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.stressorManagements.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.stressorManagementNames[i]).grid(row=i+1, column=2, sticky=W)
            
        #controls
        Button(self.page, text="Back", command=self.GetIntensityNames, width=4).grid(row=int(self.stressorManagements.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetAreaChangeNames, width=4).grid(row=int(self.stressorManagements.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-management"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetAreaChangeNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'change in habitat area\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.areaChanges.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.areaChangeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetManagementNames, width=4).grid(row=int(self.areaChanges.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetStructureChangeNames, width=4).grid(row=int(self.areaChanges.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-area"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetStructureChangeNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'change in habitat structure\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.structureChanges.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.structureChangeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetAreaChangeNames, width=4).grid(row=int(self.structureChanges.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetDisturbanceFrequencyNames, width=4).grid(row=int(self.structureChanges.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
       
        self.helpKey="criteria-structure"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetDisturbanceFrequencyNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'frequency of natural disturbance\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.disturbanceFrequencies.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.disturbanceFrequencyNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetStructureChangeNames, width=4).grid(row=int(self.disturbanceFrequencies.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetOverlapTimeNames, width=4).grid(row=int(self.disturbanceFrequencies.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-disturbance"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetOverlapTimeNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the temporal overlap levels in increasing order.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.overlapTimes.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.overlapTimeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetDisturbanceFrequencyNames, width=4).grid(row=int(self.overlapTimes.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetMortalityRateNames, width=4).grid(row=int(self.overlapTimes.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-overlap"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetMortalityRateNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'natural mortality rate\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.mortalityRates.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.mortalityRateNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetOverlapTimeNames, width=4).grid(row=int(self.mortalityRates.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecruitmentPatternNames, width=4).grid(row=int(self.mortalityRates.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-mortality"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetRecruitmentPatternNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'natural recruitment rate\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.recruitmentPatterns.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.recruitmentPatternNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetMortalityRateNames, width=4).grid(row=int(self.recruitmentPatterns.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetConnectivityNames, width=4).grid(row=int(self.recruitmentPatterns.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-recruitment"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetConnectivityNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'connectivity\' categories.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        for i in range(int(self.connectivities.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.connectivityNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetRecruitmentPatternNames, width=4).grid(row=int(self.connectivities.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecoveryTimeNames, width=4).grid(row=int(self.connectivities.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-connectivity"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetRecoveryTimeNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Categories")
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the \'age at maturity/recovery time\' categories.  Here age at maturity applies to biotic habitats and recovery time applies abiotic habitats.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6)

        for i in range(int(self.recoveryTimes.get())):
            Label(self.page, text=str(i)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.recoveryTimeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Back", command=self.GetConnectivityNames, width=4).grid(row=int(self.recoveryTimes.get())+1,column=0,sticky=W)
        Button(self.page, text="Done", command=self.EditSurvey, width=4).grid(row=int(self.recoveryTimes.get())+1,column=3,sticky=E)

        self.page.pack()
       
        self.helpKey="criteria-recovery"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()



    #data GUI
    def GetNewData(self):
        """
        
        """
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        
        if not self.importedData:
            self.habitats=StringVar()
            self.stressors=StringVar()
            self.GetHabitatStressors()
        else:
            if not self.page == None:
                self.page.pack_forget()

            self.page = Frame(self.tk)

            Label(self.page, text="Would you like to edit the data?", wraplength=self.promptWidth, justify=LEFT).grid(row=0,column=0,columnspan=4)
            Button(self.page, text="Yes", command=self.GetHabitatStressors, width=3).grid(row=1,column=1,sticky=W)
            Button(self.page, text="No", command=self.GetNewWeights, width=3).grid(row=1,column=2,sticky=E)

            self.page.pack()                  
       
        self.helpKey="data"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetHabitatStressors(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        #allow starting GUI here
        if not self.page==None:
            self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the number of habitats and stressors.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4)

        #habitats
        Label(self.page, text="Habitats").grid(row=1,column=1,sticky=E)
        Entry(self.page, textvariable=self.habitats, width=2, validate="key", validatecommand=self.vcmdInt).grid(row=1,column=2,sticky=W)

        #stressors
        Label(self.page, text="Stressors").grid(row=2,column=1,sticky=E)
        Entry(self.page, textvariable=self.stressors, width=2, validate="key", validatecommand=self.vcmdInt).grid(row=2,column=2,sticky=W)

        #controls
        Button(self.page, text="Back", command=self.GetNewData, width=4).grid(row=3,column=0,sticky=W)
        Button(self.page, text="Next",command=self.GetNewHabitatStressors, width=4).grid(row=3,column=3,sticky=E)
        self.page.pack()

        self.helpKey="data-habitats-stressors"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetNewHabitatStressors(self):
        """
        Initializes variables.
        """
        try:
            assessmentName=self.assessmentName.get()
            habitats=int(self.habitats.get())
            stressors=int(self.stressors.get())

            #change habitats if number changed or first pass and not imported data
            habitatChange=False
            if (self.habitats2 <> None and self.habitats2 <> habitats) or (self.habitats2 == None and not self.importedData):
                habitatChange=True
                habitatNames=[""]*habitats
                habitatQuality=[0]*habitats
            else:
                habitatNames=map(self.GetVarValue,self.habitatNames)
                habitatQuality=map(self.GetVarValue,self.habitatQuality)
                
            #change stressors if number changed or first pass and not imported data
            stressorChange=False
            if (self.stressors2 <> None and self.stressors2 <> stressors) or (self.stressors2 == None and not self.importedData):
                stressorChange=True
                stressorNames=[""]*stressors
                stressorQuality=[0]*stressors
                stressorIntensity=[0]*stressors
                stressorIntensityQuality=[0]*stressors
                stressorManagement=[0]*stressors
                stressorManagementQuality=[0]*stressors
                stressorDistance=[0]*stressors
            else:
                stressorNames=map(self.GetVarValue,self.stressorNames)
                stressorQuality=map(self.GetVarValue,self.stressorQuality)
                stressorIntensity=map(map(self.GetVarValue,self.stressorIntensityNames).index,map(self.GetVarValue,self.stressorIntensity))
                stressorIntensityQuality=map(self.GetVarValue,self.stressorIntensityQuality)
                stressorManagement=map(map(self.GetVarValue,self.stressorManagementNames).index,map(self.GetVarValue,self.stressorManagement))
                stressorManagementQuality=map(self.GetVarValue,self.stressorManagementQuality)
                stressorDistance=map(self.GetVarValue,self. stressorDistance)

            #change data if habitats or stressors changed or first pass and not imported data
            if habitatChange or stressorChange or ((self.habitats2 == None or self.stressors2 == None) and not self.importedData):
                areaChange=[0]*habitats*stressors
                areaChangeQuality=[0]*habitats*stressors
                structureChange=[0]*habitats*stressors
                structureChangeQuality=[0]*habitats*stressors
                disturbanceFrequency=[0]*habitats*stressors
                disturbanceFrequencyQuality=[0]*habitats*stressors
                overlapTime=[0]*habitats*stressors
                overlapTimeQuality=[0]*habitats*stressors
                mortalityRate=[0]*habitats
                mortalityRateQuality=[0]*habitats
                recruitmentPattern=[0]*habitats
                recruitmentPatternQuality=[0]*habitats
                connectivity=[0]*habitats
                connectivityQuality=[0]*habitats
                recoveryTime=[0]*habitats
                recoveryTimeQuality=[0]*habitats

                self.LoadData(assessmentName,
                              habitatNames,habitatQuality,
                              stressorNames,stressorQuality,
                              stressorIntensity,stressorIntensityQuality,
                              stressorManagement,stressorManagementQuality,
                              stressorDistance,
                              areaChange,areaChangeQuality,
                              structureChange,structureChangeQuality,
                              disturbanceFrequency,disturbanceFrequencyQuality,
                              overlapTime,overlapTimeQuality,
                              mortalityRate,mortalityRateQuality,
                              recruitmentPattern,recruitmentPatternQuality,
                              connectivity,connectivityQuality,
                              recoveryTime,recoveryTimeQuality)
                
            self.GetHabitatNames()
        except ValueError, e:
            #print e
            self.GetHabitatStressors()

    def GetHabitatNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.stressors2 = int(self.stressors.get())
        self.habitats2 = int(self.habitats.get())
        
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify names for the habitats.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=5+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Data Quality").grid(row=1,column=3, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Entry(self.page, width=20, textvariable=self.habitatNames[i], validate="key", validatecommand=self.vcmdString).grid(row=2+i,column=2)
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.habitatQuality[i]).grid(row=2+i,column=3+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=3+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetHabitatStressors, width=4).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.DuplicateHabitats, width=4).grid(row=3+int(self.habitats.get()),column=4+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-habitat-names"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def DuplicateHabitats(self):
        if len(set(map(self.GetVarValue,self.habitatNames))) < int(self.habitats.get()):
            if tkMessageBox.askyesno("Duplicate Habitats", "Are you sure you want duplicate habitat names?", default=tkMessageBox.NO, icon=tkMessageBox.WARNING):
                self.GetStressorNames()
        else:
            self.GetStressorNames()
            
    def GetStressorNames(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify names for the stressors.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=5+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Data Quality").grid(row=1,column=3, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Entry(self.page, width=20, textvariable=self.stressorNames[i]).grid(row=2+i,column=2)
            
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorQuality[i]).grid(row=2+i,column=3+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=3+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetHabitatNames, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.DuplicateStressors, width=4).grid(row=3+int(self.stressors.get()),column=4+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-stressor-names"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def DuplicateStressors(self):
        if len(set(map(self.GetVarValue,self.stressorNames))) < int(self.stressors.get()):
            if tkMessageBox.askyesno("Duplicate Stressors", "Are you sure you want duplicate stressor names?", default=tkMessageBox.NO, icon=tkMessageBox.WARNING):
                self.GetIntensity()
        else:
            self.GetIntensity()
            
    def GetIntensity(self):
        if not self.dataIndexBuilt:
            self.indexmenu.delete(2)
            self.dataIndex()
            self.weightIndex()
            
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify stressor intensity.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Intensity").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.stressorIntensity[i]]+map(lambda n: n.get(),self.stressorIntensityNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorIntensityQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetStressorNames, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetManagement, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-intensity"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetManagement(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify stressor management effectiveness.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=4+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Management").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.stressorManagement[i]]+map(lambda n: n.get(),self.stressorManagementNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorManagementQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetIntensity, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetDistance, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-management"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            

    def GetDistance(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the zone of influence for each stressor.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=5)

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Distance (m)").grid(row=1,column=3)
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            Entry(self.page, width=20, textvariable=self.stressorDistance[i]).grid(row=2+i,column=3)
                                                     
        #controls
        Button(self.page, text="Back", command=self.GetManagement, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetAreaChange, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-distance"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetAreaChange(self):
        self.habitatIndex=0
        self.GetHabitatAreaChange()

    def GetPreviousHabitatAreaChange(self):
        if self.habitatIndex == 0:
            self.GetDistance()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatAreaChange()

    def GetNextHabitatAreaChange(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatAreaChange()
        else:
            self.GetStructureChange()

    def GetHabitatAreaChangeIndex(self,n):
        self.habitatIndex=n
        self.GetHabitatAreaChange()

    def GetHabitatAreaChange(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the change in habitat area for "+self.habitatNames[self.habitatIndex].get()+" habitat when exposed to each stressor.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Change in Area").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.areaChange[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.areaChangeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.areaChangeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetPreviousHabitatAreaChange, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatAreaChange, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()
            

        self.helpKey="data-area"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            

    def GetStructureChange(self):
        self.habitatIndex=0
        self.GetHabitatStructureChange()

    def GetPreviousHabitatStructureChange(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatAreaChange()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatStructureChange()

    def GetNextHabitatStructureChange(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatStructureChange()
        else:
            self.GetDisturbanceFrequency()

    def GetHabitatStructureChangeIndex(self,n):
        self.habitatIndex=n
        self.GetHabitatStructureChange()

    def GetHabitatStructureChange(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the change in habitat structure for "+self.habitatNames[self.habitatIndex].get()+" habitat when exposed to each stressor.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Change in Structure").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.structureChange[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.structureChangeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.structureChangeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetPreviousHabitatStructureChange, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatStructureChange, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-structure"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetDisturbanceFrequency(self):
        self.habitatIndex=0
        self.GetHabitatDisturbanceFrequency()

    def GetPreviousHabitatDisturbanceFrequency(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatStructureChange()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatDisturbanceFrequency()
            
    def GetNextHabitatDisturbanceFrequency(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatDisturbanceFrequency()
        else:
            self.GetTimeOverlap()

    def GetHabitatDisturbanceFrequencyIndex(self,n):
        self.habitatIndex=n
        self.GetHabitatDisturbanceFrequency()

    def GetHabitatDisturbanceFrequency(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="For "+self.habitatNames[self.habitatIndex].get()+" habitat, specify the frequency of the natural disturbance that is most similar to each of the following stressors.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Disturbance Frequency").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.disturbanceFrequency[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.disturbanceFrequencyNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.disturbanceFrequencyQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetPreviousHabitatDisturbanceFrequency, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatDisturbanceFrequency, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-disturbance"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()


    def GetTimeOverlap(self):
        self.habitatIndex=0
        self.GetHabitatTimeOverlap()

    def GetPreviousHabitatTimeOverlap(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatDisturbanceFrequency()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatTimeOverlap()

    def GetNextHabitatTimeOverlap(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatTimeOverlap()
        else:
            self.GetMortalityRate()

    def GetHabitatTimeOverlapIndex(self,n):
        self.habitatIndex=n
        self.GetHabitatTimeOverlap()

    def GetHabitatTimeOverlap(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="For "+self.habitatNames[self.habitatIndex].get()+" habitat, specify the amount of time that the stressor and habitat come into contact.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Temporal Overlap").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.overlapTime[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.overlapTimeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.overlapTimeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetPreviousHabitatTimeOverlap, width=4).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatTimeOverlap, width=4).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-time"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetMortalityRate(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify natural mortality rate.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Mortality").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.mortalityRate[i]]+map(lambda n: n.get(),self.mortalityRateNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.mortalityRateQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetHabitatTimeOverlap, width=4).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecruitmentPattern, width=4).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-mortality"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetRecruitmentPattern(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify natural recruitment rate.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Pattern").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.recruitmentPattern[i]]+map(lambda n: n.get(),self.recruitmentPatternNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.recruitmentPatternQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

       #controls
        Button(self.page, text="Back", command=self.GetMortalityRate, width=4).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetConnectivity, width=4).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-recruitment"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetConnectivity(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify connectivity.", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Connectivity").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.connectivity[i]]+map(lambda n: n.get(),self.connectivityNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.connectivityQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetRecruitmentPattern, width=4).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecoveryTime, width=4).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-connectivity"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetRecoveryTime(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Data")
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify age at maturity (biotic habitats) or recovery time (abiotic habitats).", wraplength=self.promptWidth, justify=LEFT).grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Recovery").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.recoveryTime[i]]+map(lambda n: n.get(),self.recoveryTimeNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.recoveryTimeQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Back", command=self.GetConnectivity, width=4).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Done", command=self.GetNewWeights, width=4).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

        self.helpKey="data-recovery"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    #weights GUI
    def GetNewWeights(self):
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Weights")
        
        Label(self.page, text="Would you like to change the weighting for each criterion?", wraplength=self.promptWidth, justify=LEFT).grid(row=0,column=0,columnspan=4)
        Button(self.page, text="Yes", command=self.GetExposureWeights, width=3).grid(row=1,column=1,sticky=W)
        Button(self.page, text="No", command=self.GetNewExport, width=3).grid(row=1,column=2,sticky=E)

        self.page.pack()
       
        self.helpKey="weights"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetExposureWeights(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Weights")
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        exposure=[(10,"Spatial Overlap"),(0,"Intensity"),(1,"Management"),(5,"Temporal Overlap")]

        Label(self.page, text="Specify the weights for the following EXPOSURE criteria.", wraplength=self.promptWidth, justify=LEFT).grid(row=0,column=0,columnspan=6)

        for i,(j,name) in enumerate(exposure):
            Label(self.page, text=name).grid(row=1+i,column=1,sticky=E)
            for k in range(len(self.weightClasses)):
                Radiobutton(self.page, text="", value=k, variable=self.weights[j]).grid(row=1+i,column=3+k)

        for l,weight in enumerate(self.weightClasses):
            Label(self.page, text=str(weight), wraplength=1).grid(row=3+i,column=3+l,sticky=N)
                           
        #controls
        Button(self.page, text="Back", command=self.GetNewWeights, width=4).grid(row=4+i,column=0,sticky=W)
        Button(self.page, text="Next",command=self.GetConsequenceWeights, width=4).grid(row=4+i,column=4+k,sticky=E)
        self.page.pack()

        self.helpKey="weights-exposure"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()
            
    def GetConsequenceWeights(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Editing Weights")
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        exposure=[(2,"Area"),(3,"Structure"),(4,"Disturbance"),(6,"Mortality"),(7,"Recruitment"),(8,"Connectivity"),(9,"Recovery")]

        Label(self.page, text="Specify the weights for the following CONSEQUENCE criteria.", wraplength=self.promptWidth, justify=LEFT).grid(row=0,column=0,columnspan=6)

        for i,(j,name) in enumerate(exposure):
            Label(self.page, text=name).grid(row=1+i,column=1,sticky=E)
            for k in range(len(self.weightClasses)):
                Radiobutton(self.page, text="", value=k, variable=self.weights[j]).grid(row=1+i,column=3+k)

        for l,weight in enumerate(self.weightClasses):
            Label(self.page, text=str(weight), wraplength=1).grid(row=3+i,column=3+l,sticky=N)
                           
        #controls
        Button(self.page, text="Back", command=self.GetExposureWeights, width=4).grid(row=4+i,column=0,sticky=W)
        Button(self.page, text="Done",command=self.GetNewExport, width=4).grid(row=4+i,column=4+k,sticky=E)
        self.page.pack()

        self.helpKey="weights-consequence"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def GetNewExport(self):
        self.tk.wm_title("Marine InVEST - Habitat Risk Assessment" + " - " +self.assessmentName.get() + " - Export")
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)
        self.tk.wm_title(self.tk.wm_title()[:self.tk.wm_title().rfind("-")-1])
        
        Label(self.page, text="Would you like to export the assessment?", wraplength=self.promptWidth, justify=LEFT).grid(row=0,column=0,columnspan=4)
        Button(self.page, text="Yes", command=self.Export, width=3).grid(row=1,column=1,sticky=W)
        Button(self.page, text="No", command=self.NoExport, width=3).grid(row=1,column=2,sticky=E)

        self.page.pack()
       
        self.helpKey="export"
        if self.showHelp.get():
            if not self.helpRoot==None:
                self.helpRoot.withdraw()
            self.OnShowHelp()

    def NoExport(self):
        if tkMessageBox.askyesno("Export", "Discard changes?", default=tkMessageBox.NO, icon=tkMessageBox.WARNING):
            self.Home()


    def Export(self):
        filename=tkFileDialog.asksaveasfilename(initialfile=self.assessmentName.get()+".csv",filetypes=[("Comma Separated Values",".csv"),("All Files","*")])
        if filename=="":
            pass
        else:
            if filename.split(".")[-1] <> "csv":
                filename = filename + ".csv"
            self.WriteFile(filename)
            self.Home()

    def Import(self):
        filename=tkFileDialog.askopenfilename(filetypes=[("Comma Separated Values",".csv"),("All Files","*")])
        if filename=="":
            pass
        else:
            self.ReadFile(filename)
            self.importedData=True
            self.tk.wm_title(self.tk.wm_title() + " - " +self.assessmentName.get())
            self.habitats2=int(self.habitats.get())
            self.stressors2=int(self.stressors.get())

            #build index
            self.indexmenu.delete(0)
            self.criteriaIndex()
            self.dataIndex()
            self.weightIndex()            

            self.EditSurvey()

if __name__ == "__main__":
    HRA()
