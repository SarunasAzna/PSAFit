"""
    v0.6.1  
    This is a package  fitPmTm  which holds function used in protein melting 
    pressure-temperature phase diagram fitting program.

    Copyright (C) 2015 Baltic Institute of Advanced Technology

    Authors: Sarunas Azna, Piotras Cimmperman, Vytautas Rafanavicius

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import numpy as np
from math import exp, log, sqrt
import psaFit as bp
import qti
from PyQt4 import QtCore, QtGui
import traceback
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
import time
import sys

## Protein melting pressure-temperature phase diagram function.
#  @param dict Dictionary|Class : function arguments
#  @param dic['Pr'] Float : reference pressure value
#  @param dic['T'] Float : reference temperature value
#  @param dic['DuCp'] Float : isobaric specific heat change of protein unfolding value
#  @param dic['DuH'] Float : enthalpy change of protein unfolding value
#  @param dic['DuV'] Float : volume change of protein unfolding value
#  @param dic['DuG'] Float : Gibbs free energy change of protein unfolding value
#  @param dic['DuBeta'] Float : isothermic change of protein unfolding value
#  @param dic['DuAlpha'] Float : thermal expansion factor change of protein unfolding value
#  @param x Float : temperature value
#  @param i=0 Int : iterations parameter
#  @return Dictionary : solutions of the quadratic equation
#  @return {'Pm1'} Float : first solution value
#  @return {'Pm2'} Float : second solution value
def ellipsoidPmTm(dict, x, i=0, other_params={}):
    try:
        dic = dict.valuesdict()
    except:
        dic = dict
    try:
        Pr =                    other_params['Pr']
        T =                     other_params['T'] #this is T reference or T0

        DuCp =                  dic['DuCp']
        DuH =                   dic['DuH']
        DuG =                   dic['DuG']

        DuBeta =                dic['DuBeta']
        DuV =                   dic['DuV']
        DuAlpha =               dic['DuAlpha']


    except KeyError:
        Pr =                    other_params['Pr']
        T =                     other_params['T'] #this is T reference or T0

        DuCp =                  dic['DuCp_%i' % (i+1)]
        DuH =                   dic['DuH_%i' % (i+1)]
        DuG =                   dic['DuG_%i' % (i+1)]

        DuBeta =                dic['DuBeta_%i' % (i+1)]
        DuV =                   dic['DuV_%i' % (i+1)]
        DuAlpha =               dic['DuAlpha_%i' % (i+1)]

    a = DuBeta/2
    b = DuV + DuAlpha * (x - T)
    c = -DuCp * (x*(log(x/T) -1) + T) - (DuH - DuG)*(x/T - 1) + DuG
    D = b**2 - 4*a*c
    try:
        Pm1 = (-b - sqrt(D))/(2*a) + Pr
        Pm2 = (-b + sqrt(D))/(2*a) + Pr
    except ValueError:
        Pm1 = 0
        Pm2 = 0

    return {'Pm1': Pm1,
                    'Pm2' : Pm2,
                    'D' : D}


## This function for creating LmFit Parameters expression to Enthalpy parameter.
#  @param fitParams Lmfit Parameter Object : function parameters
#  @param x Array : x(temperature) values
#  @param yi=0 Integer : iteration parameter
#  @return fitParams Lmfit Parameter Object : function parameters
def ellipsoidEnthalpy(fitParams, x, yi=0):
    '''
    The main objective of this function is to create parameter which holds quadratic equation solutions as real numbers --> [D]iscriminant >= 0.
    Example is in http://lmfit.github.io/lmfit-py/constraints.html segment 'Using Inequality Constrains'. Delta - discriminant, y - enthalpy(DuH), x - all other parameters.
    '''
    ## Creation of discriminant parameter, [R]ight and [L]eft
    fitParams.add('DR_%i' % (yi+1),  min=0.0, vary=True)
    fitParams.add('DL_%i' % (yi+1),  min=0.0, vary=True)
    ## Creation of temperature extremum parameter.
    '''
    The right extremum point of PmLt ellipse is used  with assumpsion that left extremum will be unreached.
    '''
    fitParams.add('extrR_%i' % (yi+1),  value = float(max(x)), vary=False)
    fitParams.add('extrL_%i' % (yi+1),  value = float(min(x)), vary=False)
    aR = '(DuBeta_%i /2)' % (yi + 1)
    aL = '(DuBeta_%i /2)' % (yi + 1)
    bR = '(DuV_%i + DuAlpha_%i * (extrR_%i - T_%i ))'  % (yi + 1, yi + 1, yi + 1, yi + 1)
    bL = '(DuV_%i + DuAlpha_%i * (extrL_%i - T_%i ))'  % (yi + 1, yi + 1, yi + 1, yi + 1)
    '''
    In case of simplicity the c parameter of quadratic equation is defided like this:
    c = f - (DuH - DuG)*(x/T - 1)
    f = -DuCp * (x*(log(x/T) -1) + T) + DuG
    '''
    fR = '(-DuCp_%i * (extrR_%i *(log(extrR_%i /T_%i ) -1) + T_%i ) + DuG_%i )' % (yi + 1, yi + 1, yi + 1, yi + 1, yi + 1, yi + 1)
    fL = '(-DuCp_%i * (extrL_%i *(log(extrL_%i /T_%i ) -1) + T_%i ) + DuG_%i )' % (yi + 1, yi + 1, yi + 1, yi + 1, yi + 1, yi + 1)
    ##Expression of enthalpy is found from this equation: D = b^2 - 4ac, with [R]ight extremum
    DuHR = ' '.join(
                                    ['(', bR, '**2', '-', 'DR_%i' % (yi + 1), '-', '4', '*', aR,  '*', fR, ')', '*',
                                    'T_%i' % (yi + 1), '/','(', '4', '*', aR, '*', '(', 'T_%i' % (yi + 1),  '-',
                                    'extrR_%i' % (yi + 1), ')', ')', '+', 'DuG_%i' % (yi+1)  ]
                                    )
    fitParams.add('DuHR_%i' % (yi+1), expr=DuHR)
    ##Expression of enthalpy is found from this equation: D = b^2 - 4ac, with [L]eft extremum
    DuHL = ' '.join(
                                    ['(', bL, '**2', '-', 'DL_%i' % (yi + 1), '-', '4', '*', aL,  '*', fL, ')', '*',
                                    'T_%i' % (yi + 1), '/','(', '4', '*', aL, '*', '(', 'T_%i' % (yi + 1),  '-',
                                    'extrL_%i' % (yi + 1), ')', ')', '+', 'DuG_%i' % (yi+1)  ]
                                    )
    fitParams.add('DuHL_%i' % (yi+1), expr=DuHL)
    ##Creation of Enthalpy parameter
    '''
    There are some technical problems with two constraints:
    fitParams.add('DuH_%i' % (yi+1), expr='DuHL_%i if DR_%i > DL_%i else DuHR_%i' % (yi + 1, yi + 1, yi + 1, yi + 1))
    Fit parameter is added before fit, both (left and right) discriminants when is  equal to 0 so the program choose firs if opsion.
    '''
    fitParams.add('DuH_%i' % (yi+1), expr='DuHR_%i' % (yi +1))
    return fitParams

## Function for fitting Pm(Tm) curves (pressure-temperature phase diagram).
#  @param fitWizard Qti Table Object : template table
#  @param currentTableCols Qti Table Object : selected columns of the current table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Table Object : results params table(parameters values)
#  @return Qti Graph Object : results graph(plot curves)
def fitPmTm(tablename = ''):
    try:
        ## initializing parameters
        startTime = time.time()
        rr = [ ["DuCp", "DuH", "DuG", "DuBeta", "DuV", "DuAlpha" ],
                        ["Rsquare" ,"Pr", "T", "DuCp", "DuH", "DuG", "DuBeta", "DuV", "DuAlpha"]]
        fitvars = ["DuCp", "DuH", "DuG", "DuBeta", "DuV", "DuAlpha"]
        t = qti.app.table(tablename)
        if  not isinstance(t, qti.Table):
            t = qti.app.currentTable()
            if not isinstance(t, qti.Table):
                msgBox = QtGui.QMessageBox()
                msgBox.setText(str(e))
                ret = msgBox.exec_()
        wizard = qti.app.table("fitWizard")
        useGlobal = False # logic parameter for using globality in the module
        answer = bp.columnToArray(tablename)
        x = answer[0]
        data =  answer[1]
        fitWiz = bp.fitWizardParameters(wizard, useGlobal, fitvars)
        prams = fitWiz[0]
        other_params = fitWiz[1]
        if useGlobal:
            globalD = fitWiz[2]
        else:
            globalD = {}
        fit_params = Parameters()
        ## initiating fit parameters
        for yi, y in enumerate(data):
            for name in prams:
                fit_params.add( name + '_%i' % (yi + 1)    , value = float(prams[name].value), min = float(prams[name].min), max = float(prams[name].max), vary = prams[name].vary)

            if fit_params['DuBeta'+ '_%i' % (yi + 1)].vary == True:
                fit_params.add('eplips_%i' % (yi + 1), value = 100, min = 1e-09, vary = True)
                fit_params['DuBeta'+ '_%i' % (yi + 1)].expr = '-(eplips_%i  + DuAlpha_%i **2)*' % (yi + 1, yi + 1)+ str(other_params['T']) +'/DuCp_%i' % ( yi + 1)
            #elif fit_params['DuAlpha'+ '_%i' % (yi + 1)].vary == True:
                #fit_params.add('eplips_%i' % (yi + 1), value = 2, min = 1e-09, vary = True)
                #fit_params['DuAlpha'+ '_%i' % (yi + 1)].expr = 'sqrt(-eplips_%i  - DuBeta_%i/' % (yi + 1, yi + 1)+ str(other_params['T']) +'*DuCp_%i)' % ( yi + 1)
            #fit_params.add('deltaDuCp_%i' % (yi+1), value = -0.3218791946, max=(0-1e-09), vary=False)

            #fit_params.add('DuCp_%i' % (yi+1), expr='(deltaDuCp_%i -DuAlpha_%i **2)*T_%i /DuBeta_%i' % (yi + 1, yi + 1, yi + 1, yi + 1))
            #bp.ellipsoidEnthalpy(fit_params, x, yi)


        ## ititiating and assigning global fit parameters
        if useGlobal:
            for name in globalD:
                if globalD[name]:
                    fit_params.add( name , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
                    for yi, y in enumerate(data):
                        fit_params[name + '_%i' % (yi + 1)].expr = name
        fit_value = []
        for u  in rr[0]:
            fit_value.append(u)
            fit_value.append(fit_params[u+'_1'].value)
            fit_value.append(fit_params[u+'_1'].min)
            fit_value.append(fit_params[u+'_1'].max)
            fit_value.append(fit_params[u+'_1'].vary)
        funcname = bp.ellipsoidPmTm
        def objective(params, x, data, funcname, other_params):
            ndata, nx = data.shape
            resid = 0.0*data[:]
            for i in range(ndata):
                for k in range(nx):
                    dat = data[i, k]
                    tempres = funcname(params, x[k], i=i, other_params= other_params)
                    if k == 0 or (k > 0 and (x[k] >= x[k-1])):
                        model = tempres['Pm1']
                    elif k > 0 and (x[k] < x[k-1]):
                        model = tempres['Pm2']
                    if tempres['D'] > 0:
                        resid[i, k] = dat - model
                    else:
                        resid[i, k] = 50
            return resid.flatten()
        ## FIT
        fitStartTime = time.time()
        minimizer = minimize(objective, fit_params, args=(x, data, funcname, other_params))

        rezStartTime = time.time()

        fitRezTName = bp.fitResultElip(x, funcname, fit_params,  other_params, t, rr[0], fitTableName = "TfitResElip_1")
        a = bp.Rsquare(qti.app.table(fitRezTName))
        p_values = bp.runs_test(qti.app.table(fitRezTName))
        chi_squares = bp.chi_squares(qti.app.table(fitRezTName), \
            minimizer.nvarys)
        chi_sqr = chi_squares[0]
        red_chi = chi_squares[1]

        #for yi, i in enumerate(a):
            #fit_params.add( 'Rsquare_%i' % (yi+1)    , value = i)
        rr.append('Rsquare')
        outputTable = bp.outputTable(minimizer, t, fit_params, fit_value, rr, data, fitRezName="fitRezTName", useGlobal = useGlobal, globalD = globalD, tName = 'TfitParElip_1' )
        pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                'x_name':"<i>T<sub>m</sub></i>, K",
                'change_x': True,
                'x_scale' : 0,
                'x_num_format' : 4,
                'DotsLine' : 1,
                'x_min' : min(x),
                'x_max' : max(x),
                'graphName' : 'GfitElip_1' }
        bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

        g = qti.app.currentGraph()
        fitwiz = qti.app.table('fitWizard')
        fitwiz.showMaximized()
        fitwiz.showNormal()
        g.showMaximized()
        g.showNormal()

        ##Write rezults to log
        qti.app.updateLog(' \n ')
        qti.app.updateLog("Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + ' \n ')

        qti.app.updateLog("Input table name: " + str(t.objectName()) + ' \n ')
        qti.app.updateLog("Output table name: " + str(outputTable.objectName()) + ' \n ')
        qti.app.updateLog("Fit result table name: " + fitRezTName + ' \n ')
        qti.app.updateLog("Fit result graph name: " + g.objectName() + ' \n ')
        qti.app.updateLog("function: " + funcname.__name__ + ' \n ')
        for yi, i in enumerate(a):
            qti.app.updateLog( 'Rsquare_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(p_values):
            qti.app.updateLog( 'P_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(chi_sqr):
            qti.app.updateLog( 'Chisqr_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(red_chi):
            qti.app.updateLog( 'redChi_%i' % (yi+1) + ": " + str(i) + ' \n ')
        qti.app.updateLog(fit_report(fit_params) + ' \n ')
        endTime = time.time()
        #msgBox = QtGui.QMessageBox()
        #msgBox.setText("Done")
        qti.app.updateLog('Time: ' + str("%.4g"%( endTime - startTime) )+'s \n ' +
                                                        'Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n ' +
                                                        'Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n ' +
                                                        'Results displaying time: ' + str("%.4g"%(endTime - rezStartTime))+'s \n ' +
                                                        '--------------------------------------------------------------------------------------------------------------------  \n'
                                                        )
        #ret = msgBox.exec_()

    except Exception, e:
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Error!!!")
        if str(e) == "'NoneType' object has no attribute 'selectedColumns'" or str(e) == "integer division or modulo by zero":
            msgBox.setInformativeText("You must to select x and y columns from the table")
        else:
            msgBox.setInformativeText(str(e))
        msgBox.setDetailedText(str(traceback.format_exc()))
        ret = msgBox.exec_()
        qti.app.resultsLog().append(str(traceback.format_exc()))
