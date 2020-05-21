import math
import matplotlib.pyplot as plt
import numpy as np
import xlwt 
from xlwt import Workbook 

# mesh #size of mesh
# icl #classical inversion point
# nodes 
# ncross #number of times solution changes signs

xmax = int(input('What is the max_value of x for integration limits(x_max): '))
mesh = int(input('number of divisions you want: '))

dx = xmax/mesh #Grid Size

Excitation_state_start = int(input('From which excited state you want to calculate(Put 0 for ground, 1 for first and so on): '))
Excitation_state_end = int(input('Until which excitation state do you want to calculate: '))

def Eigen(xmax, mesh, Excitation_state_start, Excitation_state_end, dx):

    fig = plt.figure()
    pos = 1
    for state in range(Excitation_state_start,Excitation_state_end+1):

        number_of_nodes = state
        vpot = []
        x = []
        y = []
        p = []
        func = []
        X = []
        Y = []
        P = []
        V = []
        Y2 = []

        for i in range(0,mesh):
            vpot.append(0)
            x.append(0)
            y.append(0)
            p.append(0)
            func.append(0)

        ddx12 = dx*dx/12.0 #delta x square for calculations

        #Set up the potential

        for i in range(0,mesh):
            x[i] = i*dx
            vpot[i] = 0.5*x[i]*x[i]

        # Eigen Value search
        #Set the initial lower and upper bounds to the eigen value


        eup = max(vpot)
        elow = min(vpot)

        #Set the initial guess of the energy as the mid point of eup and elow

        e = 0.5*(eup + elow)
        n_iter = 1000 #Number of times to iterate so as to find the perfect energy

        for kkk in range(1,n_iter+1):

            #func is the numerov algorithm's function
            #we deteremine the position of last crossing i.e. the sign change

            func[0] = ddx12*(2.0*(vpot[0]-e))
            icl = -1
            key = 0
            for i in range(1,mesh):
                func[i] = ddx12*2.0*(vpot[i]-e)
                #Beyond the classical point of inversion this function will change sign
                if( func[i] == 0):
                    func[i] = 10**-20 ## Ensuring the value of function didn't go to zero
                if (func[i]*func[i-1] < 0 and key != 1):
                    icl = i
                    key = 1

            if (icl >= mesh-2 or icl < 1):
                for i in range(len(func)):
                    print("no classical turning point")
                    break
            else:
                for i in range(0,mesh):
                    func[i] = 1 - func[i]
            
            hnodes = number_of_nodes/2

            if(number_of_nodes%2 == 1):
                y[0] = 1
                y[1] = (0.5*(12 - 10*func[0])*y[0])/func[1]
            else:
                y[0] = 0
                y[1] = dx

            #Start Outward integration and count number of times signs change
            ncross = 0
            for i in range(1,mesh-1):
                y[i+1] = ((12 - 10*func[i])*y[i] - func[i-1]*y[i-1])/func[i+1]  #from i=0 to i = icl

                if( y[i]*y[i+1] <0):
                    ncross += 1
    
            for i in range(mesh-2,icl+1,-1):
                y[i-1] = ((12-10*func[i])*y[i] - func[i+1]*y[i+1])/func[i-1]

            if (n_iter > 1):
                if(number_of_nodes%2 == 0):
                    if (ncross > hnodes-1):  # Number of sign changes are more than number of nodes, 
                        eup = e            #hence energy is too high, lower the upper bound
                    else:
                        elow = e            # Energy lower, increase lower bound
                else:
                    if(ncross>hnodes):
                        eup=e
                    else:
                        elow = e
            
            e = 0.5 * (eup + elow) #This gives us the next trial value

            if(eup - elow < 10**-10):
                break

        y[mesh-1] = 0
        y[mesh-2] = dx

        # for i in range(mesh-2,icl+1,-1):
        #     y[i-1] = ((12-10*func[i])*y[i] - func[i+1]*y[i+1])/func[i-1] #From i=mesh-1 to i=icl+1
        #     # print(y[i-1])

        # # To compare the derivatives sign at i = icl, we need to find y[icl], y[icl-1] 
        # # from both inward and outward integration

        y_right_icl = ((12-10*func[icl+1])*y[icl+1] - func[icl+2]*y[icl+2])/func[icl]
        for i in range(icl+1,mesh):
            y[i] = y[i]*y[icl]/y_right_icl

        norm = 0
        for i in range(icl,mesh):
            p[i] = 0

        for i in range(0,icl):
            arg = (e - (x[i]**2)/2)
            if( arg > 0):
                p[i] = 1/ math.sqrt(arg)
            else:
                p[i] = 0

            norm = norm + 2*dx*p[i]

        norm = norm - dx*p[0]

        #Normalize p[x] so that Integration of p[x]dx = 1
        for i in range(0,icl):
            p[i] = p[i]/norm


        #Output arrays X , Y[x] , Y[x]**2, classical P[x], V
        power_factor=-1
        if number_of_nodes == 0:
            power_factor=1
        else:
            for i in range(0,number_of_nodes):
                power_factor *= -1 

        E = []
        for i in range(mesh-1,0,-1):
            X.append(-1*x[i])
            Y.append(power_factor*y[i] + e)
            # if(number_of_nodes%2 == 0):
            #     Y.append(-1*power_factor*y[i] + e)
            # else:
            #     Y.append(power_factor*y[i] + e)
            Y2.append(y[i]*y[i])
            P.append(p[i])
            V.append(vpot[i])
            E.append(e)

        for i in range(0,mesh):
            X.append(x[i])
            Y.append(y[i] + e)
            Y2.append(y[i]*y[i])
            P.append(p[i])
            V.append(vpot[i])
            E.append(e)

        if(state == 1):
            print("Ground State eigen value is " + str(e))
        else:
            print("Energy Eigen value of excited state " + str(state) + " is " + str(e))

        name = "WaveFunctions" + str(state) + ".xls"
        workbook = xlwt.Workbook()
        sheet = workbook.add_sheet('Sheet 1')
        for i in range(0,2*mesh-1):
            sheet.write(i,0,X[i])
            sheet.write(i,1,Y[i])
            sheet.write(i,2,Y2[i])
        workbook.save(name)

        colours = ['blue','red','green','black','yellow']
        X = np.array(X)
        Y = np.array(Y)
        # plt.plot(X,Y,linewidth=1.0, color=colours[state])
        ax = fig.add_subplot(Excitation_state_end-Excitation_state_start+1,2,pos)
        ax1 = fig.add_subplot(Excitation_state_end-Excitation_state_start+1,2,pos+1)
        pos = pos+2
        ax.plot(X,Y,linewidth=1.0, color=colours[state])
        ax.plot(X,E,linewidth=1.0, color='black', label='EigenValue i.e. Energy Levels', linestyle='dashed')
        ax.legend()
        ax.set_ylabel(r'$\psi_{%d}$'%(state))
        ax1.plot(X,Y2,linewidth=1.0, color=colours[state])
        ax1.set_ylabel(r'$\psi^2_{%d}$'%(state))
Eigen(xmax,mesh,Excitation_state_start+1, Excitation_state_end+1, dx)

plt.show()