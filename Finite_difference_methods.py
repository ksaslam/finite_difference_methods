# --------------------------------------------------------
# This codes solves the differential equation
#    (d^2)u/dx^2 + du/dx = 0
# using many different partial differential methods like
# Euler (forward and backward), Trapozoidal, Bashford, Leapfrog etc
# ---------------------------------------------------------  


import numpy as np
import matplotlib.pyplot as plt
import math as mt



final_t =8.;   # final time to solve
increment= np.array([0.001, 2*0.001, 4*0.001, 8*0.001, 16*0.001 ] )


##############################################################
# Initialization of different fields , u is displacement and v is velocity
for k in xrange(0,5):
  dt_a= increment[k];
  n= final_t * np.pi/ dt_a
  t_a= np.linspace(0,final_t* np.pi,n)
  #print len(t_a)
  v_ef= np.ones((len(t_a)))
  u_ef= np.zeros((len(t_a)))
  v_eb= np.ones((len(t_a)))
  u_eb= np.zeros((len(t_a)))
  v_trp= np.ones((len(t_a)))
  u_trp= np.zeros((len(t_a)))
  v_leapf= np.ones((len(t_a)))
  u_leapf= np.zeros((len(t_a)))
  v_bash= np.zeros((len(t_a)))
  u_bash= np.zeros((len(t_a)))
  u_back= np.zeros((len(t_a)))
  v_back= np.zeros((len(t_a)))
#########################################################################


  # Analytical solution 
  u_analy= np.zeros((len(t_a)))

  u_analy= ( 1./ mt.sqrt(2)  ) * np.sin( (mt.sqrt(2) * t_a) )
  v_analy=  np.cos( (mt.sqrt(2) * t_a) )



 # forward Euler
  v_ef[0]=1.
  u_ef[0]=0.
  for i in xrange(0,len(t_a)-1):
    matrix1= np.array([[1., dt_a], [-2.*dt_a, 1] ])
    u_ef[i+1]= matrix1[0][0]*u_ef[i]+ matrix1[0][1]* v_ef[i]
    v_ef[i+1]= matrix1[1][0]*u_ef[i]+ matrix1[1][1]* v_ef[i]


# Backward Euler method
  v_eb[0]=1.
  u_eb[0]=0.
  for i in xrange(0,len(t_a)-1):
 
   matrix1= np.array([[1., -dt_a], [2.*dt_a, 1] ]) 
   matrix1_inv= np.linalg.inv(matrix1)
   u_eb[i+1]= matrix1_inv[0][0]*u_eb[i]+ matrix1_inv[0][1]* v_eb[i]
   v_eb[i+1]= matrix1_inv[1][0]*u_eb[i]+ matrix1_inv[1][1]* v_eb[i]


#

# Trapozoidal Euler method
  v_trp[0]=1.
  u_trp[0]=0.
  for i in xrange(0,len(t_a)-1):
   matrix1= np.array([[1., -dt_a/2.], [dt_a, 1.] ]) 
   matrix2= np.array([[1., dt_a/2.], [-dt_a, 1.] ])
   matrix1_inv= np.linalg.inv(matrix1)
   vector_trp= np.array([u_trp[i],v_trp[i]])
   matrix_res = np.dot( matrix2, vector_trp)
   u_trp[i+1]= matrix1_inv[0][0]*matrix_res[0]+ matrix1_inv[0][1]* matrix_res[1]
   v_trp[i+1]= matrix1_inv[1][0]*matrix_res[0]+ matrix1_inv[1][1]* matrix_res[1]

# bashford method
  u_bash[0]= 0.
  v_bash[0]= 1.
  #u_bash[1]= 0.
  #v_bash[1]= 1.
  u_bash[-1] = u_bash[0]- dt_a * v_bash[0]
  v_bash[-1]= v_bash[0]- 2*dt_a*u_bash[0] 
  for i in xrange(0,len(t_a)-1):
    u_bash[i+1]=  dt_a*( 3./2. * v_bash[i]- 1./2. *v_bash[i-1]  )+ u_bash[i]
    v_bash[i+1]= dt_a*( 3./2. * -2.* u_bash[i]- 1./2. *-2. *u_bash[i-1]  )+ v_bash[i]

  
# Backward differences
  matrix_3= np.array([[3.,-2.*dt_a], [4.*dt_a ,3.] ])
  matrix3_inv_back= np.linalg.inv(matrix_3)
  u_back[0]= 0.
  v_back[0]= 1.
  u_back[-1] = u_back[0]- dt_a * v_back[0]
  v_back[-1]= v_back[0]- 2*dt_a*u_back[0] 
  for i in xrange(0,len(t_a)-1):
    v_backector1= np.array( [4*u_back[i]- u_back[i-1], 4*v_back[i]- v_back[i-1] ] )
    u_back[i+1]= matrix3_inv_back[0][0]*v_backector1[0]+ matrix3_inv_back[0][1]* v_backector1[1]
    v_back[i+1]= matrix3_inv_back[1][0]*v_backector1[0]+ matrix3_inv_back[1][1]* v_backector1[1]



#leapfrog method
  v_leapf[0]=1.
  u_leapf[0]=0.
  for i in xrange(0,len(t_a)-1):
	 matrix1= np.array([[1., dt_a], [-2.*dt_a, 1.-2.*(dt_a)**2] ])
	 u_leapf[i+1]= matrix1[0][0]*u_leapf[i]+ matrix1[0][1]* v_leapf[i]
	 v_leapf[i+1]= matrix1[1][0]*u_leapf[i]+ matrix1[1][1]* v_leapf[i]

   

  # This part below is for plotting purpose only ################################# 

  if (k==0):
   t_a_dt1= t_a
   u_analy_dt1= u_analy
   v_analy_dt1= v_analy
   u_eb_dt1= u_eb
   v_eb_dt1= v_eb
   u_ef_dt1= u_ef
   v_ef_dt1= v_ef
   u_trp_dt1=u_trp
   v_trp_dt1= v_trp
   u_leapf_dt1= u_leapf
   v_leapf_dt1= v_leapf
   u_bash_dt1= u_bash
   v_bash_dt1= v_bash
   u_back_dt1= u_back
   v_back_dt1= v_back


  if (k==1):
   t_a_dt2=t_a
   u_analy_dt2= u_analy
   u_eb_dt2= u_eb
   v_eb_dt2= v_eb
   u_ef_dt2= u_ef
   v_ef_dt2= v_ef
   u_trp_dt2=u_trp
   v_trp_dt2= v_trp
   u_leapf_dt2= u_leapf
   v_leapf_dt2= v_leapf  
   u_bash_dt2= u_bash
   v_bash_dt2= v_bash   
   u_back_dt2= u_back
   v_back_dt2= v_back

  if (k==2):
   t_a_dt3= t_a
   u_analy_dt3= u_analy
   u_eb_dt3= u_eb
   v_eb_dt3= v_eb
   u_ef_dt3= u_ef
   v_ef_dt3= v_ef
   u_trp_dt3=u_trp
   v_trp_dt3= v_trp
   u_leapf_dt3= u_leapf
   v_leapf_dt3= v_leapf
   u_bash_dt3= u_bash
   v_bash_dt3= v_bash
   u_back_dt3= u_back
   v_back_dt3= v_back  
 
  if (k==3):
   t_a_dt4= t_a 
   u_analy_dt4= u_analy
   u_eb_dt4= u_eb
   v_eb_dt4= v_eb
   u_ef_dt4= u_ef
   v_ef_dt4= v_ef
   u_trp_dt4=u_trp
   v_trp_dt4= v_trp
   u_leapf_dt4= u_leapf
   v_leapf_dt4= v_leapf 
   u_bash_dt4= u_bash
   v_bash_dt4= v_bash
   u_back_dt4= u_back
   v_back_dt4= v_back

  if (k==4):
   t_a_dt5= t_a 
   u_analy_dt5= u_analy
   u_eb_dt5= u_eb
   v_eb_dt5= v_eb
   u_ef_dt5= u_ef
   v_ef_dt5= v_ef
   u_trp_dt5=u_trp
   v_trp_dt5= v_trp
   u_leapf_dt5= u_leapf
   v_leapf_dt5= v_leapf 
   u_bash_dt5= u_bash
   v_bash_dt5= v_bash
   u_back_dt5= u_back
   v_back_dt5= v_back

 
################## Simple Plots #######################################

plt.figure(1)
plt.plot(t_a_dt1, u_ef_dt1, label= 'forward Euler') 
plt.plot( t_a, u_analy, label='Analytical')
plt.title('forward Euler method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)   


plt.figure(2)
plt.plot(t_a_dt1, u_eb_dt1, label= 'Euler Backward')
plt.plot( t_a, u_analy, label='Analytical')   
plt.title('Backward Euler method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)  

plt.figure(3)
plt.plot(t_a_dt1, u_trp_dt1, label= 'Trapoziodal')
plt.plot( t_a, u_analy, label='Analytical')   
plt.title('trapozoidal method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)  

plt.figure(4)
plt.plot(t_a_dt1, u_leapf_dt1, label= 'leapfrog')
plt.plot( t_a, u_analy, label='Analytical')   
plt.title('leapfrog method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)  

plt.figure(5)
plt.plot(t_a_dt1, u_bash_dt1, label='Adams-Bashforth') 
plt.plot( t_a, u_analy, label='Analytical')
plt.title('Adam- bashford method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)  
#plt.show()
### 5 and 6 plots are missing 

plt.figure(6)
plt.plot(t_a_dt1, u_back_dt1, label= 'backward difference') 
plt.plot( t_a, u_analy, label='Analytical')
plt.title('Backward Difference method - dt = 0.001')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right')
plt.ylim(-1.2, 1.2)  
#plt.show()

############################ Phase plots #######################################

plt.figure(7)
plt.plot(u_ef_dt1, v_ef_dt1,label='forward Euler')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical')   
plt.title(' phase pane plot- forward Euler method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)
#plt.show()

plt.figure(8)
plt.plot(u_eb_dt1, v_eb_dt1, label='Backward Euler')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical')   
plt.title(' phase pane plot- Backward Euler method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)

#plt.show()

plt.figure(9)
plt.plot(u_trp_dt1, v_trp_dt1, label= 'Trapoziodal')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical')
plt.title(' phase pane plot- Trapoziodal method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)

#plt.show()

plt.figure(10)
plt.plot(u_leapf_dt1, v_leapf_dt1, label='leapfrog')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical') 
plt.title(' phase pane plot-  Leapfrog method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)

#plt.show()

plt.figure(11)
plt.plot(u_bash_dt1, v_bash_dt1,label= 'adams-bashforth')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical') 
plt.title(' phase pane plot-  Adam- Bashforth method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)

plt.figure(12)
plt.plot(u_back_dt1, v_back_dt1, label= 'backward difference')
plt.plot( u_analy_dt1, v_analy_dt1 , label='Analytical') 
plt.title(' phase pane plot-  Backward Difference method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right')
plt.ylim(-1.5, 1.5)
plt.xlim(-1,1)







### Calculating the error #####################

 # error values for backward Euler
er_eb1= abs(u_analy_dt1 - u_eb_dt1)
max_er_eb1= max(er_eb1)# finding the max values 
er_eb2= abs(u_analy_dt2- u_eb_dt2)
max_er_eb2= max(er_eb2)
er_eb3= abs(u_analy_dt3- u_eb_dt3)
max_er_eb3= max(er_eb3)
er_eb4= abs(u_analy_dt4- u_eb_dt4)
max_er_eb4= max(er_eb4)
er_eb5= abs(u_analy_dt5- u_eb_dt5)
max_er_eb5= max(er_eb5)
max_eb= np.array([max_er_eb1,max_er_eb2,max_er_eb3,max_er_eb4, max_er_eb5])

er_ef1= abs(u_analy_dt1 - u_ef_dt1)
max_er_ef1= max(er_ef1)# finding the max values 
er_ef2= abs(u_analy_dt2- u_ef_dt2)
max_er_ef2= max(er_ef2)
er_ef3= abs(u_analy_dt3- u_ef_dt3)
max_er_ef3= max(er_ef3)
er_ef4= abs(u_analy_dt4- u_ef_dt4)
max_er_ef4= max(er_ef4)
er_ef5= abs(u_analy_dt5- u_ef_dt5)
max_er_ef5= max(er_ef5)
max_ef= np.array([max_er_ef1,max_er_ef2,max_er_ef3,max_er_ef4, max_er_ef5])


er_bash1= abs(u_analy_dt1 - u_bash_dt1)
max_er_bash1= max(er_bash1)# finding the max values 
er_bash2= abs(u_analy_dt2- u_bash_dt2)
max_er_bash2= max(er_bash2)
er_bash3= abs(u_analy_dt3- u_bash_dt3)
max_er_bash3= max(er_bash3)
er_bash4= abs(u_analy_dt4- u_bash_dt4)
max_er_bash4= max(er_bash4)
er_bash5= abs(u_analy_dt5- u_bash_dt5)
max_er_bash5= max(er_bash5)
max_bash= np.array([max_er_bash1,max_er_bash2,max_er_bash3,max_er_bash4, max_er_bash5])

er_leapf1= abs(u_analy_dt1 - u_leapf_dt1)
max_er_leapf1= max(er_leapf1)# finding the max values 
er_leapf2= abs(u_analy_dt2- u_leapf_dt2)
max_er_leapf2= max(er_leapf2)
er_leapf3= abs(u_analy_dt3- u_leapf_dt3)
max_er_leapf3= max(er_leapf3)
er_leapf4= abs(u_analy_dt4- u_leapf_dt4)
max_er_leapf4= max(er_leapf4)
er_leapf5= abs(u_analy_dt5- u_leapf_dt5)
max_er_leapf5= max(er_leapf5)
max_leapf= np.array([max_er_leapf1,max_er_leapf2,max_er_leapf3,max_er_leapf4, max_er_leapf5])

er_trp1= abs(u_analy_dt1 - u_trp_dt1)
max_er_trp1= max(er_trp1)# finding the max values 
er_trp2= abs(u_analy_dt2- u_trp_dt2)
max_er_trp2= max(er_trp2)
er_trp3= abs(u_analy_dt3- u_trp_dt3)
max_er_trp3= max(er_trp3)
er_trp4= abs(u_analy_dt4- u_trp_dt4)
max_er_trp4= max(er_trp4)
er_trp5= abs(u_analy_dt5- u_trp_dt5)
max_er_trp5= max(er_trp5)
max_trp= np.array([max_er_trp1,max_er_trp2,max_er_trp3,max_er_trp4, max_er_trp5])

er_back1= abs(u_analy_dt1 - u_back_dt1)
max_er_back1= max(er_back1)# finding the max values 
er_back2= abs(u_analy_dt2- u_back_dt2)
max_er_back2= max(er_back2)
er_back3= abs(u_analy_dt3- u_back_dt3)
max_er_back3= max(er_back3)
er_back4= abs(u_analy_dt4- u_back_dt4)
max_er_back4= max(er_back4)
er_back5= abs(u_analy_dt5- u_back_dt5)
max_er_back5= max(er_back5)
max_back= np.array([max_er_back1,max_er_back2,max_er_back3,max_er_back4, max_er_back5])






############################### phase plot for all the methods together ############

plt.figure(15)
plt.plot( u_trp_dt1, v_trp_dt1, u_trp_dt2 , v_trp_dt2, u_trp_dt3, v_trp_dt3, u_trp_dt4, v_trp_dt4,u_trp_dt5, v_trp_dt5,u_analy_dt1,v_analy_dt1 )  
plt.title(' phase pane plot-  Trapoidal method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.ylim(-1.1, 1.1)
plt.xlim(-1,1)
#plt.show()


plt.figure(16)
plt.plot( u_ef_dt1, v_ef_dt1, label= 'dt= 0.001')
plt.plot(u_ef_dt2 , v_ef_dt2, u_ef_dt3, v_ef_dt3, u_ef_dt4, v_ef_dt4,u_ef_dt5, v_ef_dt5,u_analy_dt1,v_analy_dt1 )  
plt.title(' phase pane plot-  Forward Euler method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right',prop={'size':12})
plt.ylim(-1.5, 1.5)
plt.xlim(-1.2,1.2)
plt.show()

plt.figure(17)
plt.plot( u_eb_dt1, v_eb_dt1, label= 'dt=0.001')
plt.plot(u_eb_dt2 , v_eb_dt2, u_eb_dt3, v_eb_dt3, u_eb_dt4, v_eb_dt4,u_eb_dt5, v_eb_dt5,u_analy_dt1,v_analy_dt1 ,u_analy_dt1,v_analy_dt1)  
plt.title(' phase pane plot-  Backward Euler method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)')
plt.legend(loc='upper right',prop={'size':12})
plt.ylim(-1.1, 1.1)
plt.xlim(-.8,.8)

plt.figure(18)
plt.plot( u_leapf_dt1, v_leapf_dt1, u_leapf_dt2 , v_leapf_dt2, u_leapf_dt3, v_leapf_dt3, u_leapf_dt4, v_leapf_dt4,u_leapf_dt5, v_leapf_dt5,u_analy_dt1,v_analy_dt1 ) 
plt.title(' phase pane plot-  Leapfrog method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)') 

#plt.ylim(-1.1, 1.1)
#plt.xlim(-.8,.8)

plt.figure(18)
plt.plot( u_bash_dt1, v_bash_dt1, u_bash_dt2 , v_bash_dt2, u_bash_dt3, v_bash_dt3, u_bash_dt4, v_bash_dt4,u_bash_dt5, v_bash_dt5,u_analy_dt1,v_analy_dt1 ) 
plt.title(' phase pane plot-  Adam Bashford method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)') 

plt.ylim(-1.1, 1.1)
plt.xlim(-.8,.8)

plt.figure(37)
plt.plot( u_back_dt1, v_back_dt1, u_back_dt2 , v_back_dt2, u_back_dt3, v_back_dt3, u_back_dt4, v_back_dt4,u_back_dt5, v_back_dt5,u_analy_dt1,v_analy_dt1 ) 
plt.title(' phase pane plot-  Backward Difference method ')
plt.xlabel('u(t)')
plt.ylabel('v(t)') 

plt.ylim(-1.1, 1.1)
plt.xlim(-.8,.8)

#######################################################All plots of displacement 

plt.figure(19)
plt.plot(t_a_dt1, u_trp_dt1, label= 'dt=0.001')
plt.plot(t_a_dt2, u_trp_dt2 , label= '2*dt')
plt.plot(t_a_dt3, u_trp_dt3, label= '4*dt')
plt.plot(t_a_dt4, u_trp_dt4, label= '8*dt')
plt.plot(t_a_dt5,u_trp_dt5, label= '16*dt')
plt.plot(t_a_dt1,u_analy_dt1 , label= 'Analytical')
plt.title('Trapozoidal method ')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right',prop={'size':12})
plt.xlim(0,31)
plt.show()

plt.figure(20)
plt.plot(t_a_dt1, u_eb_dt1, label= 'dt=0.001')
plt.plot( t_a_dt2, u_eb_dt2,label= '2*dt')
plt.plot(t_a_dt3, u_eb_dt3,label= '4*dt')
plt.plot(t_a_dt4, u_eb_dt4, label= '8*dt')
plt.plot( t_a_dt5,u_eb_dt5,label= '16*dt' )   
plt.plot(t_a_dt1,u_analy_dt1 ,label= 'Analytical')
plt.title('Backward Euler method ')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend(loc='upper right',prop={'size':12})
plt.xlim(0,31)
plt.show()

plt.figure(21)
plt.plot(t_a_dt1, u_ef_dt1, label= 'dt=0.001')
plt.plot(t_a_dt2, u_ef_dt2 ,label= '2*dt')
plt.plot(t_a_dt3, u_ef_dt3,label= '4*dt')
plt.plot(t_a_dt4, u_ef_dt4,label= '8*dt')
plt.plot(t_a_dt5,u_ef_dt5, label= '16*dt')   
plt.title('Forward Euler method')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.xlim(0,31)
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(22)
plt.plot(t_a_dt1, u_leapf_dt1, label= 'dt=0.001')
plt.plot(t_a_dt2, u_leapf_dt2 ,label= '2*dt')
plt.plot(t_a_dt3, u_leapf_dt3,label= '4*dt')
plt.plot(t_a_dt4, u_leapf_dt4, label= '8*dt')
plt.plot(t_a_dt5,u_leapf_dt5, label= '16*dt')   
plt.title('Leapfrog method ')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.plot(t_a_dt1,u_analy_dt1 , label= 'Analytical')
plt.xlim(0,31)
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(23)
plt.plot(t_a_dt1, u_bash_dt1, label= 'dt= 0.001')
plt.plot(t_a_dt2, u_bash_dt2 , label= '2*dt')
plt.plot(t_a_dt3, u_bash_dt3,label= '4*dt')
plt.plot(t_a_dt4, u_bash_dt4, label= '8*dt')
plt.plot(t_a_dt5,u_bash_dt5, label= '16*dt')
plt.title('Adam- Bashford method ')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.plot(t_a_dt1,u_analy_dt1 , label= 'Analytical')
plt.xlim(0,31)
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(24)
plt.plot(t_a_dt1, u_back_dt1,label= 'dt= 0.001')
plt.plot(t_a_dt2, u_back_dt2 , label= '2*dt')
plt.plot(t_a_dt3, u_back_dt3, label= '4*dt')
plt.plot(t_a_dt4, u_back_dt4, label= '8*dt')
plt.plot(t_a_dt5,u_back_dt5, label= '16*dt')  
plt.title('Backward Difference method ')
plt.xlabel('time')
plt.ylabel('u(t)')
plt.plot(t_a_dt1,u_analy_dt1 , label= 'Analytical')
plt.xlim(0,31)
plt.legend(loc='upper right',prop={'size':12})
plt.show()

############################## Error plots of all the methods #################################

plt.figure(38)

plt.semilogy(t_a_dt1, er_eb1, label='backward Euler')
plt.semilogy(t_a_dt1, er_ef1, label='Forward Euler')
plt.semilogy(t_a_dt1, er_bash1, label='Adam Bashford')
plt.semilogy(t_a_dt1, er_back1, label='backward Difference')
plt.semilogy(t_a_dt1, er_trp1, label='Trapozoidal')
plt.semilogy(t_a_dt1, er_leapf1, label='Leapfrog')
plt.legend(loc=2,prop={'size':8})
plt.legend(loc='upper left')
#t_a_dt1, er_ef1, t_a_dt1,er_leapf1, t_a_dt1, er_trp1, t_a_dt1, er_bash1,t_a_dt1,er_back1)
plt.title(' Error plot-  All methods ')
plt.xlabel('time')
plt.ylabel('Error')
plt.ylim(1e-9, 1)
plt.show()

#############################################################################################

plt.figure(25)
plt.semilogy(t_a_dt1, er_eb1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_eb2, t_a_dt3,er_eb3, t_a_dt4, er_eb4, t_a_dt5, er_eb5)
plt.title(' Error plot-  Backward Euler method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(26)
plt.semilogy(t_a_dt1, er_ef1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_ef2, t_a_dt3,er_ef3, t_a_dt4, er_ef4, t_a_dt5, er_ef5)
plt.title(' Error plot-  Forward Euler method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(27)
plt.semilogy(t_a_dt1, er_bash1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_bash2, t_a_dt3,er_bash3, t_a_dt4, er_bash4, t_a_dt5, er_bash5)
plt.title(' Error plot-  Adams-Bashford method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(28)
plt.semilogy(t_a_dt1, er_trp1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_trp2, t_a_dt3,er_trp3, t_a_dt4, er_trp4, t_a_dt5, er_trp5)
plt.title(' Error plot-  Trapozoidal method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()

plt.figure(29)
plt.semilogy(t_a_dt1, er_leapf1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_leapf2, t_a_dt3,er_leapf3, t_a_dt4, er_leapf4, t_a_dt5, er_leapf5)
plt.title(' Error plot-  Leapfrog method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()


plt.figure(30)
plt.semilogy(t_a_dt1, er_back1, label= 'dt=0.001')
plt.semilogy(t_a_dt2, er_back2, t_a_dt3,er_back3, t_a_dt4, er_back4, t_a_dt5, er_back5)
plt.title(' Error plot-  Backward Difference method ')
plt.xlabel('time')
plt.ylabel('Error')
plt.legend(loc='upper right',prop={'size':12})
plt.show()

######################################## Convergence plots



# Convergence plot
plt.figure(31)
plt.loglog(increment, max_eb, label= ' backward Euler')
plt.loglog(increment, max_ef, label= 'forward Euler')
plt.plot(increment, max_trp, label=' trapozoidal')
plt.plot(increment, max_leapf, label= ' leapfrog')
plt.plot(increment, max_bash, label='adams-bashforth')
plt.loglog(increment, max_back, label='backward difference')
plt.title(' Convergence plot-  All method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
plt.legend(loc=2,prop={'size':12})
#plt.legend(loc='upper left')
plt.ylim(1e-3, 1)
plt.xlim(1e-3, 2e-2)
plt.show()

plt.figure(32)
plt.plot(increment, max_ef)
plt.title(' Convergence plot-  Forward Euler method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
#plt.show()

plt.figure(33)
plt.plot(increment, max_trp)
plt.title(' Convergence plot- Trapozoidal method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
#plt.show()

plt.figure(34)
plt.plot(increment, max_leapf)
plt.title(' Convergence plot-  Leapfrog method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
#plt.show()

plt.figure(35)
plt.plot(increment, max_bash)
plt.title(' Convergence plot-   Adams-Bashford method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
#plt.show()

plt.figure(36)
plt.loglog(increment, max_back)
plt.title(' Convergence plot-  Backward Difference method ')
plt.xlabel('delta t')
plt.ylabel('max. Error')
plt.show()


