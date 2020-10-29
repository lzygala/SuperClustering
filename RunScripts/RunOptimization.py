import os, sys, subprocess

runAnalyzer = True
runOptimization = True
runCaloShapePlotter = True

dataFile = 'root://eoscms.cern.ch///eos/cms/store/user/lzygala/ReReco_OptPFRecHitThresholds/4Gamma_sample_OptPFRHT.root'       #input for EnvelopeAnalyzer
analyzerFile = 'Output/PlotFiles/EnvelopeAnalyzer_Output.root'   #output from EnvelopeAnalyzer & input for EnvelopeOptimizer
optimizerFile = 'Output/PlotFiles/EnvelopeOptimizer_Output.root'  #output from EnvelopeOptimizer & input for CaloShapePlotter

localRegionParameterFile = 'Output/ParameterLists/localParameters.txt'
averagedParameterFile = 'Output/ParameterLists/averagedParameters.txt'
dPhiParameterFile = 'Output/ParameterLists/dPhiParameters.txt'

logEBins = 15
minLogE = -2
maxLogE = 1

seedEtaBins = 15
minSeedEta = 0      #abs
minSeedEta = 3      #abs

caloClusterShapeDEtaDistBins = 120
dPhiWindowEtaBins = 4
dPhiWindowETDistBins = 200

#optimizer options
error = 2
ringRejection = 1
binRejection = 1

def argConversion(arg):
    if type(arg) in [float, int]:
        return str(arg)
    elif type(arg) in [str]:
        return '"' + arg + '"'
    else:
        return '"' + str(arg) + '"'

def argConnection(arglist=[]):
    if not (arglist is [] or arglist is None):
        shellargv = list(arglist) 
        arglist_conv = list()
        for arg in arglist:
            arglist_conv.append(argConversion(arg))
        return '('+','.join(arglist_conv)+')'
    else:
        return ''

def runMacro(macroName, arglist=None, splash=False, interprete=False, batch=True):
    shellCommand = ['root']
    if interprete is False:
        shellCommand.append("-q")
    if splash is False:
        shellCommand.append("-l")
    if batch is True:
        shellCommand.append("-b")
    shellCommand.append(macroName+argConnection(arglist))
    print("Run Macro", shellCommand)
    a = subprocess.Popen(shellCommand)
    return a

def main():
    analyzerArgList = [dataFile, #input file 
                       analyzerFile] #output file


    optimizerArgList = [analyzerFile, #input file,
                        optimizerFile, #output file
                        localRegionParameterFile,
                        averagedParameterFile,
                        dPhiParameterFile,
                        error,
                        ringRejection,
                        binRejection]

    plotterArgList = [optimizerFile, #input file
                      localRegionParameterFile,
                      averagedParameterFile,
                      dPhiParameterFile] 

    if runAnalyzer:
        print("Running Analyzer")
        analyzer = runMacro('Analyzers/EnvelopeAnalyzer.C++', analyzerArgList)
        analyzer.wait()

    if runOptimization:
        print("Running Optimizer")
        optimization = runMacro('Analyzers/EnvelopeOptimizer.C++', optimizerArgList)
        optimization.wait()

    if runCaloShapePlotter:
        print("Running Plotter")
        plotter = runMacro('Plotters/CaloShapePlotter.C++', plotterArgList)
        plotter.wait()

if __name__ == "__main__":
    main()