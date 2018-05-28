#!/usr/bin/python

# compute vertical levels using 2 tanh algorithm
import numpy as np
from math import *
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider, Button, RadioButtons


#functions
def findcoef(zacr,zacr2,zkth,zkth2 ) :
    alpha1 = zacr  * log ( cosh( (1. -zkth  ) / zacr  ) )
    alpha2 = zacr2 * log ( cosh( (1. -zkth2 ) / zacr2 ) )
    alpha3 = zacr  * log ( cosh( (jpk-zkth  ) / zacr  ) )
    alpha4 = zacr2 * log ( cosh( (jpk-zkth2 ) / zacr2 ) )
    beta1  =         tanh(       (1.5-zkth  ) / zacr  )
    beta2  =         tanh(       (1.5-zkth2 ) / zacr2 )
    beta3  =         tanh(       (jpk-0.5-zkth  ) / zacr  )
    beta4  =         tanh(       (jpk-0.5-zkth2 ) / zacr2 )
    # declare matrix
    A=np.matrix( [[1.,1.,alpha1, alpha2], [1., jpk, alpha3, alpha4 ], [ 0., 1., beta1, beta2],[ 0.,1, beta3, beta4]] )
    B=np.matrix( [[h0],[hmax],[dzmin],[dzmax] ] )
#   [zsur, za0, za1, za2]=np.linalg.solve(A,B)
    return np.linalg.solve(A,B)

def printcoef() : #zacr,zkth,zsur,za0,za1,za2,zkth2,zacr2) :
    print "ppsur     =",  zsur
    print "ppaa0     = ", za0
    print "ppa1      = ", za1
    print "ppkth     =",  zkth
    print "ppacr     = ", zacr
    print "ldbleranh = ", ldbletan
    print "ppa2      = ", za2
    print "ppkth2    =",  zkth2
    print "ppacr2    = ", zacr2

def printall(val) : 
    printcoef()
    jk=0
    while ( jk < jpk ) :
        print "%4d %9.2f %9.2f %9.3f %9.3f " % (jk+1, gdepw_t[jk], gdept_t[jk],e3w_t[jk],e3t_t[jk])
        jk=jk+1
   

def depfunc( zsur, za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan ) :
   zw = np.arange(1  ,jpk+1  )
   zt = np.arange(1.5,jpk+1.5)
   ztmpw1 =  zsur + za0 * zw + za1 * zacr * np.log ( np.cosh( (zw-zkth ) / zacr  ) )
   if ( ldbletan ) :
      ztmpw2 = za2 * zacr2* np.log ( np.cosh( (zw-zkth2) / zacr2 ) )
      ztmpw = ztmpw1 + ztmpw2

   ztmpt1 =  zsur + za0 * zt + za1 * zacr * np.log ( np.cosh( (zt-zkth ) / zacr  ) )
   if ( ldbletan  ) :
      ztmpt2 = za2 * zacr2* np.log ( np.cosh( (zt-zkth2) / zacr2 ) )
      ztmpt = ztmpt1 + ztmpt2

   return ( ztmpw1, ztmpw2, ztmpw, ztmpt1, ztmpt2, ztmpt )

def e3func( za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan ) :
   zw = np.arange(1  ,jpk+1  )
   zt = np.arange(1.5,jpk+1.5)
   ztmpw1 = za0            + za1        * np.tanh(       (zw-zkth ) / zacr  )
   if ( ldbletan ) :
      ztmpw2  =              za2        * np.tanh(       (zw-zkth2) / zacr2 )
      ztmpw   = ztmpw1 + ztmpw2 

   ztmpt1 = za0            + za1        * np.tanh(       (zt-zkth ) / zacr  )
   if ( ldbletan ) :
      ztmpt2  =              za2        * np.tanh(       (zt-zkth2) / zacr2 )
      ztmpt   = ztmpt1 + ztmpt2 

   return ( ztmpw1, ztmpw2, ztmpw, ztmpt1, ztmpt2, ztmpt )

def update(val):
   global gdepw_t, gdept_t, e3w_t, e3t_t
   global zkth, zacr, zkth2, zacr2, zsur, za0, za1, za2
   zkth = skth.val
   zacr = sacr.val
   zkth2 = skth2.val
   zacr2 = sacr2.val
   solution=findcoef(zacr,zacr2,zkth,zkth2)
   zsur =  solution [0,0]
   za0  =  solution [1,0]
   za1  =  solution [2,0]
   za2  =  solution [3,0]

#   printcoef() #zacr,zkth,zsur,za0,za1,za2,zkth2,zacr2) 
   gdepw_0, gdepw_2, gdepw_t, gdept_0, gdept_2, gdept_t = depfunc(  zsur, za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan )
   e3w_0,   e3w_2,   e3w_t,   e3t_0,   e3t_2,   e3t_t   = e3func(         za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan )

   l1.set_xdata(gdepw_t)
   l1.set_ydata(e3t_t)
   l2.set_ydata(gdepw_t)
   l3.set_ydata(e3t_t)
   l4.set_ydata(e3t_0)
   l5.set_ydata(e3t_2)

   fig.canvas.draw_idle()
#  fig1.canvas.draw_idle()
#  fig2.canvas.draw_idle()
#  fig3.canvas.draw_idle()
#  fig4.canvas.draw_idle()

def reset(event):
    skth.reset()
    sacr.reset()
    skth2.reset()
    sacr2.reset()
    update

def finish(event):
    print "Ciao"
    quit()
   
#====================================================================================
h0    = 0.
hmax  = 6000.
dzmin = 1.
dzmax = 50
jpk   = 300
ldbletan = True
# coef to be set before (default value)
zacr  = 38.4998412469  #   80 # 7.0*3
zkth  = 317.238187329  # 184.2121644 #15.35101370000000*jpk/75.*3
zacr2 = 121.356963487 #  13.0*3
zkth2 = 31.5541059316 # 48.029893720000*jpk/75.

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.05, bottom=0.25)
axcolor = 'lightgoldenrodyellow'

axkth   = plt.axes([0.25, 0.16, 0.65, 0.02], axisbg=axcolor)
axacr   = plt.axes([0.25, 0.12, 0.65, 0.02], axisbg=axcolor)
axkth2  = plt.axes([0.25, 0.08, 0.65, 0.02], axisbg=axcolor)
axacr2  = plt.axes([0.25, 0.04, 0.65, 0.02], axisbg=axcolor)

skth2   = Slider(axkth2, 'kth2', 0.1, 400.0, valinit=zkth2)
skth    = Slider(axkth , 'kth ', 0.1, 400.0, valinit=zkth )
sacr    = Slider(axacr,  'acr',  0.1, 150,   valinit=zacr )
sacr2   = Slider(axacr2, 'acr2', 0.1, 150.0, valinit=zacr2)

A=np.ones( (4,4) )
B=np.ones( (4,1 ) ) 

solution=findcoef(zacr,zacr2,zkth,zkth2)
zsur =  solution [0,0]
za0  =  solution [1,0]
za1  =  solution [2,0]
za2  =  solution [3,0]

gdepw_0, gdepw_2, gdepw_t, gdept_0, gdept_2, gdept_t = depfunc(  zsur, za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan )
e3w_0,   e3w_2,   e3w_t,   e3t_0,   e3t_2,   e3t_t   = e3func(         za0, za1, zkth, zacr, za2, zkth2, zacr2 , ldbletan )

printcoef() #zacr,zkth,zsur,za0,za1,za2,zkth2,zacr2) 

fig1=plt.subplot(2,2,1)
l1,=plt.plot(gdepw_t, e3t_t, 'g^')
gr1=plt.grid(True)
plt.axis([0, 6000, 0, 60])

fig2=plt.subplot(2,2,2)
l2, =plt.plot(gdepw_t ,'r^')
gr2 =plt.grid(True)

fig3=plt.subplot(2,2,3)
l3,=plt.plot(e3t_t ,'b^')
gr3=plt.grid(True)

fig4=plt.subplot( 2,2,4 )
l4,l5 =plt.plot(e3t_0 ,'b^', e3t_2, 'b+' )
gr4=plt.grid(True)

resetax = plt.axes([0.02, 0.075, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

resetax2 = plt.axes([0.02, 0.025, 0.1, 0.04])
button2 = Button(resetax2, 'Quit', color='red', hovercolor='0.575')

resetax3 = plt.axes([0.02, 0.125, 0.1, 0.04])
button3 = Button(resetax3, 'Print', color='blue', hovercolor='0.575')

button.on_clicked(reset)
button2.on_clicked(finish)
button3.on_clicked(printall)

skth2.on_changed(update)
skth.on_changed(update)
sacr.on_changed(update)
sacr2.on_changed(update)

plt.show()

