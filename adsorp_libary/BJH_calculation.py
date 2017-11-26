import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .utilis import *
 
#------------------------- class ---------------------------------------

class BJH_method():
    def __init__(self, pore_size_min = 17, gas_type='N2', ):
        self.pore_size_min = pore_size_min
        self.gas_type = gas_type
        #self.p_res, self.q_res = se
            

    def fit_transform(self, p, q):
        #calculate porosity 
        self.p = p
        self.q = q
        self.vpore_total, self.vpore_micro, self.vpore_meso = _get_porosity(self.p,self.q,gas_type=self.gas_type)
        Davg, LP, Dp, dV_desorp, k = BJH_main(self.p, self.q, self.gas_type)
        
        # only select pore size >= pore_size_min
        id_valid_pore = Davg>self.pore_size_min
        self.Davg, self.LP, self.Dp, self.dV_desorp  = Davg[id_valid_pore],\
            LP[id_valid_pore], Dp[id_valid_pore], dV_desorp[id_valid_pore]
        self.k = len(self.Davg)-1
        print(self.k)

        self.Vp, self.Vp_ccum, self.Vp_dlogD = _get_volume(self.Davg, self.LP, self.Dp)
		
        return self
		
    def plot_isotherm(self):
        figure = plt.figure()
        legend = []
        legend_raw, = plt.plot(self.p_raw,self.q_raw,'ko-',label='raw iso')
        legend.append(legend_raw)
        if self.use_pressure:
            legend_fix, = plt.plot(self.p, self.q, 'r.',label='fixed iso')
            legend.append(legend_fix)
        plt.legend(handles=legend,loc=4)
        plt.grid()

    def plot_BJH_psd(self,plot_type = 'incremental'):
        figure = plt.figure()
        legend = []
        if plot_type ==   'incremental':
            plt.title('PSD by incremental pore volume')
            legend_incre, = plt.semilogx(self.Davg[1:], self.Vp[1:], 'go-', label='incremental')
            legend.append(legend_incre)
        plt.legend(handles=legend)
        plt.grid()

				
def _get_volume(Davg,LP,Dp):
    """ """

    Vp = np.pi*LP*(Davg/2.0)**2 *10**(-16) # return Vp vector[cm^3/g]
    Vp_ccum = np.add.accumulate(Vp)
    Vp_dlogD = np.zeros(len(Vp)-1)
    for i in range(len(Vp)-1):
        Vp_dlogD[i] = Vp[i]/ np.log10(Dp[i]/Dp[i+1])
    return Vp,Vp_ccum,Vp_dlogD
	
def _get_porosity(p, q, gas_type='N2'):
    '''Calculate the porosity of total, meso and micro

    Args:


    Return: 
        v_pore_total,
        v_pore_micro,
        v_pore_meso
    '''

    p_res, q_res = restrict_isotherm(p, q, Pmin=0.3, Pmax=0.999)
    p0,q0 = p_res[0],q_res[0]
    p1,q1 = p_res[-1],q_res[-1]
    print(p0,q0,p1,q1)
    gas_const = get_gas_constant(gas_type)
    vpore_total = q1*gas_const['Vmol'] / 22414.0
    vpore_micro = q0*gas_const['Vmol'] / 22414.0
    vpore_meso = vpore_total - vpore_micro
    return vpore_total,vpore_micro,vpore_meso

	
#---------------- main calculation function

def BJH_main(p,Q,gas_type='N2'):
    """Calculate the pore size distribution from given pressure and adsorption quantity,

    Notice: All the pressure here are relative pressure

    Terms: 
        p_rels: np.array, relative pressure, 
        VL: liquid equivalent volume, unit [cm^3/g]
        Rc: Kelvin radius, unit [A] 
        gas_const['A']: adsorbate property factor
        gas_const['Vmol']: the liquid molar volume, unit [cm^3/mol]
        Davg: np.array, Weighted avarage pore diameter, unit[A]
        del_tw: float, change in thickness of the wall due to desorption from previously opened pores
        LP: np.array, length of previously opened pores, 
        Vd: np.array, total volume of gas desorbed from walls of previously opened pores
        dV_desorp: np.array, total volume of gas desorped at current pressure interval
        SAW: float, total surface of walls exposed so far ( come from previously opened pores), unit [cm^2/g]
        Vc: np.array, volume desorped from newly opened pores, unit [cm^3/g]
        Pavg: np.array, relative pressure corresponding to Davg; no unit, range [0,1]
        tw_avg: np.array, thickness of the adsorbed layer at pressure of Pavg 
        del_td: float, the decrease in thickness of the wall layer by desorption from the walls of new pores 
            during decrease from Pavg to end of current pressure interval 
        CSA_c: np.array, cross-sectional area of the newly opened pores, unit [cm^2/g]
        Dp: np.array, diameters corresponding to the ends of the intervals, unit[A]
        i_step: int, interval number, i=1, from P1 to P2 (P1 is the highest pressure point)
        j: int, index for each previous interval, range [1,k+1)
        k: int, total number of pressure intervals in which new pores have been found. 

    Args:
        p: numpy.array, a list of relative pressure, 
        Q: numpy.array, a list of quantity adsorption, unit [mol/g] 
        gas_type: str, the type of gas, choice from ['N2','Ar']


    Returns:
        Davg, np.array, a list of average pore diameter
        LP,
        Dp,
        dV_desorp,
        k
    """
    gas_const = get_gas_constant(gas_type)
    
    # insert the pressure=0, adsorption = 0 point 
    p = insert_zero(p)
    Q = insert_zero(Q)
    # make the isotherm in reverse order
    p_reverse = p[::-1]
    Q_reverse = Q[::-1]
    p_rels = np.zeros(len(p_reverse))
    q_ads  = np.zeros(len(p_reverse))
    p_rels[:] = p_reverse
    q_ads[:] = Q_reverse
    #print('old p_rels',p_rels,q_ads)
    p_rels,q_ads

    # Convert adsorption quantity into liquid equivalent volume  
    VL = q_ads*gas_const['Vmol'] / 22414.0 # 22414 cm^3 STP

    n_point = len(p_rels)
    n_step = n_point-1
    Vd = np.zeros(n_step)
    Vc = np.zeros(n_step)
    dV_desorp= np.zeros(n_step)
	
	# Status: the adsorption status for each step, it is initiated with 0, 
	# If it is 1, then no new pores are created. If it is 2, then addtional new 
	# pores are created. 
    status = np.zeros(n_step) 
    tw = np.zeros(n_point)
    #print('old tw/status',n_point,n_step,tw,status)

    # not using first index
    # from p_rels to Lp, all have length of number of initial points + 2. 
    p_rels, q_ads, tw, =insert_zero(p_rels), insert_zero(q_ads), insert_zero(tw)
    VL,Vd, Vc, dV_desorp, status = insert_zero(VL),insert_zero(Vd),insert_zero(Vc),insert_zero(dV_desorp),insert_zero(status)
    #print('new p_rels,q_ads',p_rels,q_ads)
    #print('VL',VL)
    #print('new tw/status',tw,status)
    # define other vector 
    Rc, Davg, Pavg= np.zeros(len(Vd)), np.zeros(len(Vd)), np.zeros(len(Vd))
    tw_avg, CSA_c, LP = np.zeros(len(Vd)), np.zeros(len(Vd)), np.zeros(len(Vd))
    #print('tw_avg',tw_avg)
    # end of parameter preparation

    # initiation of calculation
    Rc[1]  = kelvin_radius(p_rels[1],gas_const)
    tw[1] = thickness_Harkins_Jura(p_rels[1])
    #print('Rc[1]/tw[1]',Rc[1],tw[1])
    k=0
    for istep in range(1,n_step): # modified from n_step +1 to n_step
        #print('\nistep/nstep',istep,n_step)
        status[istep]= 0 
        # Calculate the thickness of next pressure step(istep+1)
        if istep == n_step: # don't have this
            tw[istep+1]=0
        else:
            tw[istep+1] = thickness_Harkins_Jura(p_rels[istep+1])

        # a) determine Vd 
        del_tw = tw[istep] - tw[istep+1] # change of thickness
        #print('del_tw',del_tw)
        #print('Vd',Vd)
        Vd[istep] = _get_CSA_a(del_tw,Davg,LP,k,istep,n_step)
        #print('Vd vs Vd_test',Vd[istep],Vd_test[istep])

        # b) check Vd with true desorption
        dV_desorp[istep] = VL[istep] - VL[istep+1]
        #print('dV_desorp',dV_desorp[istep])
        if Vd[istep] >= dV_desorp[istep]: 
        # True case 1: Vd is larger than current increment of volume desorbed dV_desorp[istep], 
        # desorption from walls only is occuring
            status[istep] = 1
            #print('too large check case ',status[istep])
            #print('too large dV_desorp ',dV_desorp[istep])
            SAW = 0
            for j in range(1,k+1):
                SAW += np.pi*LP[j]*Davg[j] * 10**(-8)
            del_tw = dV_desorp[istep]/SAW  * 10**(8) # simplified version
            #print('SAW,new del_tw',SAW,del_tw)
        else:
        # case 2: Vd < dV_desorp[istep], addtional desorption comes from new pores
            status[istep] = 2 # case 2: normal case
            #print('normal check case ',status[istep])
            Vc[istep] = dV_desorp[istep]- Vd[istep]
            #print('dV_desorp,Vc',dV_desorp[istep],Vc[istep])
            k += 1 # total number of intervals that create new pores + 1
            #print('n_pore',k)
            Rc[k+1]  = kelvin_radius(p_rels[k+1],gas_const)
            Davg[k] = 2* (Rc[k]+Rc[k+1]) *Rc[k]*Rc[k+1] / (Rc[k]**2+Rc[k+1]**2) # mathmatical average
            Pavg[k] = np.exp(-2*gas_const['A'] / Davg[k])
            tw_avg[k] = thickness_Harkins_Jura(Pavg[k])
            del_td = tw_avg[k] - tw[istep+1]
            CSA_c[k] = np.pi*(Davg[k]/2.0+del_td)**2 *10**(-16)
            LP[k] = Vc[istep]/CSA_c[k]
            #print('Rc',Rc[k],Rc[k+1])
            #print('Vc,Davg,Pavg',Vc[istep],Davg[k],Pavg[k],tw_avg[k],CSA_c[k],LP[k])

        # c) updated pore diameter to the end pressure
        if status[istep]==2: # case 2 update current new pore diameter, as we have new pore created
            #print('updated new pore')
            Davg[k] += 2*del_td
        # no matter new pore created or not, updated previous diameter
        for j in range(1,k):
            #print('updated old pore',1,k-1)
            Davg[j] += 2*del_tw
        for j in range(1,k+1):
            Rc[j] += del_tw
        #print('Davg,Rc,LP',Davg,Rc,LP)

        # for test
        #print('istep',istep)
        #print('Rc',Rc)
        #temp1 = Davg.dot(LP)
        #Vp = np.pi*LP*(Davg/2.0)**2 *10**(-16)
        #Vp_cum = sum(Vp)
        #desorp_cum = sum(dV_desorp)
        #print('sum of Davg*LP*PI',temp1)
        #print('Vp_cum,total_desorp',Vp_cum,desorp_cum)

    #print(Davg)
    Dp  = 2*Rc
    return Davg,LP,Dp,dV_desorp,k
	
def _get_CSA_a(del_tw,Davg,LP,k,istep,n_step):
    # if it is the first step, no previous pore created
    if k==0 and istep < n_step: 
        Vd_istep = 0
    # if it is last step, no new pore will be created
    elif istep == n_step: 
        Vd_istep = 9999
    # calculate Vd 
    else: 
        #print('determine Vd >> 3 has old pore')
        Vd_istep =0
        CSA_a=np.zeros(k)
        CSA_a = insert_zero(CSA_a)
        for j in range(1,k+1):
            #CSA_a[j] = np.pi*((Rc[j]+ del_tw)**2-Rc[j]**2) *10**(-16)
            CSA_a[j] = np.pi*((Davg[j]/2.0+ del_tw)**2-(Davg[j]/2.0)**2) *10**(-16) # this one works better
            Vd_istep += LP[j]*CSA_a[j]
    return Vd_istep