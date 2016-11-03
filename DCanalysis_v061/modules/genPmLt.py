"""
    v0.6.1  
    This is a package  genPmLt  which holds function used in protein dose 
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

## Function for calculating derivative characteristics.
#  @param x 1D Numpy Array : list of x(Ligand concentration) Float values
#  @param data 2D Numpy Array : list of y(melting pressure) Float values
#  @return dic Dictionary : derivative characteristics
#  @return dic['hMax_*'] Float : maximum derivative value
#  @return dic['xHmax_*'] Float : x at maximum derivative value
#  @return dic['sHalf_*'] Float : half width at half maximum
#  @return dic['hEnd_*'] Float : derivative value at the end of the curve
#  @return dic['MhMax_*'] String : success message of hMax
#  @return dic['MxHmax_*'] String : success message of xHmax
#  @return dic['MsHalf_*'] String : success message of sHalf
#  @return dic['MhEnd_*'] String : success message of hEnd
def derivPar(x, data):
    ## messages of parameters calculations success
    answer = bp.smartData(data)
    ## editing data curves - values with false brent method are eliminated
    sdata = answer[0]
    sx = x[:answer[1]]
    ## calculating logarithmic derivative of
    derivative  = bp.logDerivative(sx, sdata)
    answer = bp.smartData(data)
    dic = {}
    for yi, y in enumerate(derivative):
        dic.update({'MhMax_%i' % (yi + 1):'OK'})
        dic.update({'MxHmax_%i' % (yi + 1):'OK'})
        dic.update({'MsHalf_%i' % (yi + 1):'OK'})
        dic.update({'MhEnd_%i' % (yi +1 ):'OK'})
        dic.update({'hMax_%i' % (yi + 1): max(y)})
        indexMax = y.tolist().index(dic['hMax_%i' % (yi + 1)])
        dic.update({'xHmax_%i' % (yi + 1) : sx[indexMax]})
        array = y[:indexMax]
        hHalf = bp.findNearest(array, dic['hMax_%i' % (yi+1)]/2)
        indexHalf = y.tolist().index(hHalf)
        xHhalf = sx[indexHalf]
        dic.update({'sHalf_%i' % (yi + 1) : (dic['xHmax_%i' % (yi + 1)] - xHhalf)})
        dic.update({'hEnd_%i' % (yi + 1) : y[-1]})
        if dic['hMax_%i' % (yi+1)] == y[-1]:
            message = "The maximum of derivative value might be unreached"
            if y[-1] == y[-2]:
                message = "There is a possible linear growth in the end of the curve. "
            dic['MhMax_%i' % (yi + 1)] = message
            dic['MxHmax_%i' % (yi + 1)] = message
        if (indexMax - indexHalf) < 6:
            dif = indexMax - indexHalf
            dic['MsHalf_%i' % (yi + 1)] = "There is just %i samples between maximum and half of the maximum derivative value. "  \
                                    "It is recommended to look at the graphical output. " %dif
        if len(sx) != len(x):
            dic['MhEnd_%i' % (yi + 1)] = "Curve has been generated from " + "{:.2e}".format(sx[0]) + " to " + "{:.2e}".format(sx[-1]) \
                                     + " (instead of " + "{:.2e}".format(x[-1]) + ") of x values."
    return dic




## Function for showing logarithmic derivative parameters in QtiPlot table ([N]ative and [U]nfolded protein binding to ligand model).
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : table with derivative parameters
def derivatParamsNU():
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam = qti.app.table("SimParam")
    output_table_name = "TderParNU_1"
    log_inform = {
        'Action' : 'NU model curves differentiation parameters obtainment',
        'Initial Parameters table' : t.objectName(),
        'Simulation Parameters table' : SimParam.objectName()}
    #                                                               Code. Do not change following text
    rr = bp.colNames(SimParam)
    d = bp.initParamsDict(t)
    x = bp.xGen(d)
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        if r == False:
            break
        else:
            for i in r:
                param = i
                ## generation of data curves
                data = bp.PmLtFullSim(x, d, param=param, SimParam=SimParam, other_params = d)
                answer = bp.derivPar(x, data)
                bp.outTabDeriv(answer, param = param, initParams = d, simParams = qti.app.table("SimParam"), tName = "TderParNU_1" )
                log_inform['Output table'] = output_table_name
                bp.results_logger(log_inform)
## Function for showing logarithmic derivative parameters in QtiPlot table ([N]ative protein binding to ligand model).
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : table with derivative parameters
def derivatParamsN():
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam = qti.app.table("SimParam")
    output_table_name = "TderParN_1"
    log_inform = {
        'Action' : 'N model curves differentiation parameters obtainment',
        'Initial Parameters table' : t.objectName(),
        'Simulation Parameters table' : SimParam.objectName()}
    # Code. Do not change following text
    rr = bp.colNames(SimParam)
    
    ## erase ligand binding to protein unfolded state parameters from options
    for i in rr:
        if i[1:3] == 'bu':
            rr.remove(i)
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        ## erase ligand binding to protein unfolded state parameters from initial parameters
        for i in ['DbuV', 'Kbu', 'DbuBeta']:
            d.pop(i, None)
        x = bp.xGen(d)
        if r == False:
            break
        else:
            for i in r:
                param = i
                ## generation of data curves
                data = bp.PmLtUpper(x, d, param=param, SimParam=SimParam, other_params=d)
                answer = bp.derivPar(x, data)
                bp.outTabDeriv(answer, param=param, initParams=d, simParams=qti.app.table("SimParam"), tName=output_table_name)
                log_inform['Output table'] = output_table_name
                bp.results_logger(log_inform)



## Function for generating Pm(Lt) curves ([N]ative and [U]nfolded protein binding to ligand model).
#  @param noiseAmplitude Float : amplitude of white noise
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def generPmLtNU(noiseAmplitude):
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam=qti.app.table("SimParam")
    log_inform = {
        'Action' : 'NU model curves generation',
        'Initial Parameters table' : t.objectName(),
        'Simulation Parameters table' : SimParam.objectName(),
        'Noise Amplitude': str(noiseAmplitude)+ " MPa"}
    pdc ={
    'y_name' : "<i>P<sub>m</sub></i>, MPa",
        'x_name' : "<i>L<sub>t</sub></i>, M",
        'change_x': True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'Dots_on' : False}
    #                                                               Code. Do not change following text
    rr = bp.colNames(SimParam)
    d = bp.initParamsDict(t)
    x = bp.xGen(d)
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        if r == False:
            break
        else:
            for i in r:
                param = i
                log_inform['Simulated Parameter'] = param
                data = bp.PmLtFullSim(x, d, param=param, SimParam=SimParam, other_params = d)
                answer = bp.smartData(data)
                sdata = answer[0]
                # add noise
                sdata = bp.addNoise(sdata, noiseAmplitude)
                sx = x[:answer[1]]
                pdc.update({'x_min' : sx[0],
                        'x_max': sx[-1],
                        'graphName': 'G'+param+'NU_1'})
                output_table = bp.generatorResultTable(sx, sdata, d, param=param, outTableName = ('T'+param+'NU_1'), SimParam=SimParam)
                output_graph = bp.Geditor(dict=pdc, fromTable=True, table=output_table)
                log_inform.update({
                    'Output Table' : output_table.objectName(),
                    'Output Graph' : output_graph.objectName()
                    })
                bp.results_logger(log_inform)


## Function for generating Pm(Lt) curves ([N]ative protein binding to ligand model).
#  @param noiseAmplitude Float : amplitude of white noise
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def generPmLtN(noiseAmplitude):
    #                                Parameters initialization
    t = qti.app.table("initParams")
    SimParam = qti.app.table("SimParam")
    log_inform = {
        'Action' : 'N model curves generation',
        'Initial Parameters table' : t.objectName(),
        'Simulation Parameters table' : SimParam.objectName(),
        'Noise Amplitude': str(noiseAmplitude)+ " MPa"}
    pdc ={
        'y_name' : "<i>P<sub>m</sub></i>, MPa",
        'x_name' : "<i>L<sub>t</sub></i>, M",
        'change_x' : True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'Dots_on' : False}

    #                                                               Code. Do not change following text
    rr = bp.colNames(SimParam)
    ## erase ligand binding to protein unfolded state parameters from options
    for i in rr:
        if i[1:3] == 'bu':
            rr.remove(i)
    ## erase ligand binding to protein unfolded state parameters from initial parameters
    for k in range(100):
        r = bp.makeButtons(rr)
        d = bp.initParamsDict(t)
        ## erase ligand binding to protein unfolded state parameters from initial parameters
        for i in ['DbuV', 'Kbu', 'DbuBeta']:
            d.pop(i, None)
        x = bp.xGen(d) 
        if r == False:
            break
        else:
            for i in r:
                param = i
                log_inform['Simulated Parameter'] = param
                data = bp.PmLtUpper(x, d, param=param, SimParam=SimParam, other_params=d)
                answer = bp.smartData(data)
                sdata = answer[0]
                # add noise
                sdata = bp.addNoise(sdata, noiseAmplitude)
                sx = x[:answer[1]]
                pdc.update({'x_min' : sx[0],
                        'x_max': sx[-1],
                        'graphName': 'G'+param+'N_1'})
                output_table = bp.generatorResultTable(sx, sdata, d, param=param, SimParam=SimParam, outTableName = 'T'+param+'N_1')
                output_graph = bp.Geditor(dict=pdc, fromTable=True, table=output_table)
                log_inform.update({
                    'Output Table' : output_table.objectName(),
                    'Output Graph' : output_graph.objectName()
                    })
                bp.results_logger(log_inform)
       
            



## Function for calculating logarithmic derivative ([N]ative and [U]nfolded protein binding to ligand model).
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def logDerivatNU():
    input_table = qti.app.currentTable()
    output_table_name = "TderivNU_1"
    a = bp.logDerivTable(tName=output_table_name , table=input_table)
    pdc ={ 
        'y_name':"<i>P<sub>m</sub></i>, MPa",
        'x_name':"<i>L<sub>t</sub></i>, M",
        'y2_name':"d<i>P<sub>m</sub></i> / d[log<sub>10</sub><i>(L<sub>t</sub>)</i>]",
        'change_x': True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'DotsLine' : 0,
        'Dots_on' : False,
        'x_max' : 100000,
        'twoYAxis' : True,
        'changeVerticalAxis' : 2,
        'graphName': 'GderivNU_1'}
    bp.Geditor(dict=pdc, fromTable=True, table=a)
    log_inform = {
        'Action' : 'NU model curves differentiation',
        'Input Table' : input_table.objectName(),
        'Output Table' : output_table_name,
        'graphName': 'GderivNU_1'
    }
    bp.Geditor(dict=pdc, fromTable=True, table=a)
    bp.results_logger(log_inform)

## Function for calculating logarithmic derivative ([N]ative protein binding to ligand model).
#  @param SimParam Qti Table Object : template table
#  @param initParams Qti Table Object : template table
#  @return Qti Table Object : results table(data curves)
#  @return Qti Graph Object : results graph(plot curves)
def logDerivatN():
    input_table = qti.app.currentTable()
    output_table_name = "TderivN_1"
    a = bp.logDerivTable(tName=output_table_name , table=input_table)
    pdc ={
        'y_name':"<i>P<sub>m</sub></i>, MPa",
        'x_name':"<i>L<sub>t</sub></i>, M",
        'y2_name':"d<i>P<sub>m</sub></i> / d[log<sub>10</sub><i>(L<sub>t</sub>)</i>]",
        'change_x': True,
        'x_scale' : 1,
        'x_num_format' : 5,
        'DotsLine' : 0,
        'Dots_on' : False,
        'x_max' : 100000,
        'twoYAxis' : True,
        'changeVerticalAxis' : 2,
        'graphName' : 'GderivN_1'}
    log_inform = {
        'Action' : 'N model curves differentiation',
        'Input table' : input_table.objectName(),
        'Output table' : output_table_name,
        'Output graph': 'GderivN_1'
    }
    bp.Geditor(dict=pdc, fromTable=True, table=a)
    bp.results_logger(log_inform)
