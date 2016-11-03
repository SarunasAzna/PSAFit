"""
    v0.6.1  
    This is a package  fitPmLt  which holds function used in protein dose
    curves fitting programs.

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

## Function for fitting Pm(Lt) curves ([N]ative protein binding to ligand model).
#  @param fitWizard Qti Table Object : template table
#  @param currentTableCols Qti Table Object : selected columns of the current table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Table Object : results params table(parameters values)
#  @return Qti Graph Object : results graph(plot curves)
def fitPmLtN(tablename):
    try:
        ## initializing parameters
        startTime = time.time()
        rr = [ ["DuG", "DuV", "DuBeta", "DbnV","DbnBeta", "Mt", "Kbn"],
                        ["Rsquare" ,"DuG", "DuV", "DuBeta", "DbnV","DbnBeta", "Mt", "Kbn", "xa", "xb", "errTol", "max_iter", "xtol", "R", "T", "Pr", "Kbu", "DbuV", "DbuBeta"]]
        fitvars = ["DuG", "DuV", "DuBeta", "Kbn", "DbnV","DbnBeta", "Mt"]
        t = qti.app.table(tablename)
        if  not isinstance(t, qti.Table):
            t = qti.app.currentTable()
            if not isinstance(t, qti.Table):
                msgBox = QtGui.QMessageBox()
                msgBox.setText(str(e))
                ret = msgBox.exec_()
        wizard = qti.app.table("fitWizard")
        useGlobal = True # logic parameter for using globality in the module
        answer = bp.columnToArray(tablename)
        x = answer[0]
        data =  answer[1]
        if data.shape[0]>14:
            rows = data.shape[0]
        else:
            rows = 15
        fitWiz = bp.fitWizardParameters(wizard, useGlobal, fitvars)
        prams = fitWiz[0]
        other_params = fitWiz[1]
        if useGlobal:
            globalD = fitWiz[2]
        else:
            globalD = {}
        fit_params = Parameters()
        ## initiating fit parameters
        for name in prams:
            for yi, y in enumerate(data):
                fit_params.add( name + '_%i' % (yi + 1)    , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
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
        funcname = bp.LtUpper
        def objective(params, x, data, funcname, other_params):
            ndata, nx = data.shape
            resid = 0.0*data[:]
            for i in range(ndata):

                for k in range(nx):
                    resid[i, k] = data[i, k] - bp.brent(funcname, x[k], params, i, other_params)
            return resid.flatten()
        ## FIT
        fitStartTime = time.time()
        minimizer = minimize(objective, fit_params, args=(x, data, funcname, other_params))
        rezStartTime = time.time()
        fitRezTName = bp.fitResult(x, funcname, fit_params, other_params, t, rr[0], brent = True, fitTableName = "TfitResN_1")

        a = bp.Rsquare(qti.app.table(fitRezTName))
        p_values = bp.runs_test(qti.app.table(fitRezTName))
        chi_squares = bp.chi_squares(qti.app.table(fitRezTName), \
            minimizer.nvarys)
        chi_sqr = chi_squares[0]
        red_chi = chi_squares[1]

        #for yi, i in enumerate(a):
                #fit_params.add( 'Rsquare_%i' % (yi+1)    , value = i)
        rr.append('Rsquare')
        outputTable = bp.outputTable(minimizer, t, fit_params, fit_value, rr, data, fitRezName=fitRezTName, useGlobal = useGlobal, globalD = globalD, tName = 'TfitParamN_1' )

        pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                'x_name':"<i>L<sub>t</sub></i>, M",
                'change_x': True,
                'x_min' : x[0],
                'x_max' : x[-1],
                'x_scale' : 1,
                'x_num_format' : 5,
                'DotsLine' : 1,
                'graphName' : 'GfitN_1' }
        bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

        g = qti.app.currentGraph()
        fitwiz = qti.app.table('fitWizard')
        fitwiz.showMaximized()
        fitwiz.showNormal()
        g.showMaximized()
        g.showNormal()
        ##Write rezults to log
        qti.app.updateLog('\n')
        qti.app.updateLog("Date: " +time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
        qti.app.updateLog("Input table name: " + str(t.objectName()) + '\n')
        qti.app.updateLog("Output table name: " + str(outputTable.objectName()) + '\n')
        qti.app.updateLog("Fit result table name: " + fitRezTName + '\n')
        qti.app.updateLog("Fit result graph name: " + g.objectName() + '\n')
        qti.app.updateLog("function: " + funcname.__name__ + '\n')
        for yi, i in enumerate(a):
            qti.app.updateLog( 'Rsquare_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(p_values):
            qti.app.updateLog( 'P_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(chi_sqr):
            qti.app.updateLog( 'Chisqr_%i' % (yi+1) + ": " + str(i) + ' \n ')
        for yi, i in enumerate(red_chi):
            qti.app.updateLog( 'redChi_%i' % (yi+1) + ": " + str(i) + ' \n ')
        qti.app.updateLog(fit_report(fit_params) + '\n')
        endTime = time.time()
        qti.app.updateLog('Time: ' + str("%.4g"%( endTime - startTime) )+'s \n ' +
                                                        'Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n' +
                                                        'Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n' +
                                                        'Results displaying time: ' + str("%.4g"%(endTime - rezStartTime))+'s \n' +
                                                        '------------------------------------------------------------------------ \n'
                                                        )

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



## Function for fitting Pm(Lt) curves ([N]ative and [U]nfolded protein binding to ligand model).
#  @param fitWizard Qti Table Object : template table
#  @param currentTableCols Qti Table Object : selected columns of the current table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Table Object : results params table(parameters values)
#  @return Qti Graph Object : results graph(plot curves)
def fitPmLtNU(tablename):
    try:
        ## initializing parameters
        startTime = time.time()
        rr = [ ["DuG", "DuV", "DuBeta", "Kbn","DbnV","DbnBeta", "Kbu","DbuV","DbuBeta","Mt"],
                        ["Rsquare" ,"DuG", "DuV", "DuBeta", "Kbn","DbnV","DbnBeta", "Kbu","DbuV","DbuBeta","Mt", "xa", "xb", "errTol", "max_iter", "xtol", "R", "T", "Pr"]]
        fitvars = ["DuG", "DuV", "DuBeta", "Kbn","DbnV","DbnBeta", "Kbu","DbuV","DbuBeta","Mt"]
        t = qti.app.table(tablename)
        if  not isinstance(t, qti.Table):
            t = qti.app.currentTable()
            if not isinstance(t, qti.Table):
                msgBox = QtGui.QMessageBox()
                msgBox.setText(str(e))
                ret = msgBox.exec_()
        wizard = qti.app.table("fitWizard")
        useGlobal = True # logic parameter for using globality in the module
        answer = bp.columnToArray(tablename)
        x = answer[0]
        data =  answer[1]

        if data.shape[0]>14:
            rows = data.shape[0]
        else:
            rows = 15
        fitWiz = bp.fitWizardParameters(wizard, useGlobal, fitvars)
        prams = fitWiz[0]
        other_params = fitWiz[1]
        if useGlobal:
            globalD = fitWiz[2]
        else:
            globalD = {}
        fit_params = Parameters()
        ## initiating fit parameters
        for name in prams:
            for yi, y in enumerate(data):
                fit_params.add( name + '_%i' % (yi + 1)    , value = prams[name].value, min = prams[name].min, max = prams[name].max, vary = prams[name].vary)
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

        funcname = bp.PmLtFull

        def objective(params, x, data, funcname, other_params):
            ndata, nx = data.shape
            resid = 0.0*data[:]
            for i in range(ndata):
                for k in range(nx):
                    resid[i, k] = data[i, k] - funcname(x[k], params, i, other_params)
            return resid.flatten()
        ## FIT
        fitStartTime = time.time()
        minimizer = minimize(objective, fit_params, args=(x, data, funcname, other_params))
        rezStartTime = time.time()
        funcnameR = bp.LtFull
        fitRezTName = bp.fitResult(x, funcnameR, fit_params, other_params, t, rr[0], brent = True, fitTableName = "TfitResNU_1")

        a = bp.Rsquare(qti.app.table(fitRezTName))
        p_values = bp.runs_test(qti.app.table(fitRezTName))
        chi_squares = bp.chi_squares(qti.app.table(fitRezTName), \
            minimizer.nvarys)
        chi_sqr = chi_squares[0]
        red_chi = chi_squares[1]
        #for yi, i in enumerate(a):
                #fit_params.add( 'Rsquare_%i' % (yi+1)    , value = i)
        rr.append('Rsquare')
        outputTable = bp.outputTable(minimizer, t, fit_params, fit_value, rr, data, fitRezName = fitRezTName, useGlobal = useGlobal, globalD = globalD, tName = 'TfitParamNU_1' )
        pdc ={ 'y_name':"<i>P<sub>m</sub></i>, MPa",
                'x_name':"<i>L<sub>t</sub></i>, M",
                'change_x': True,
                'x_min' : x[0],
                'x_max' : x[-1],
                'x_scale' : 1,
                'x_num_format' : 5,
                'DotsLine' : 1,
                'graphName' : 'GfitNU_1' }
        bp.Geditor(dict=pdc, fromTable=True, table=qti.app.table(fitRezTName))

        g = qti.app.currentGraph()
        fitwiz = qti.app.table('fitWizard')
        fitwiz.showMaximized()
        fitwiz.showNormal()
        g.showMaximized()
        g.showNormal()
        ## Write rezults to log
        qti.app.updateLog(' \n ')
        qti.app.updateLog("Date: " + time.strftime("%Y-%m-%d %H:%M:%S")  + ' \n ')
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
        qti.app.updateLog(fit_report(fit_params) + '\n')
        endTime = time.time()
        qti.app.updateLog('Time: ' + str("%.4g"%( endTime - startTime) )+'s \n ' +
                                                        'Fit time: ' + str(  "%.4g"%(rezStartTime - fitStartTime))+'s \n' +
                                                        'Data initiating time: ' + str("%.4g"%(fitStartTime - startTime))+'s \n' +
                                                        'Results displaying time: ' + str("%.4g"%(endTime - rezStartTime))+'s \n'
                                                        )

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
