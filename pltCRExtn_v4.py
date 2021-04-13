import numpy as np
import matplotlib.pyplot as plt

f0 = np.loadtxt("./f_aermix.txt")

fig = plt.figure(figsize=(6.1,5))
ax=fig.add_subplot(111)
#plt.scatter(filedata[:,1], filedata[:,2], c='navy', alpha=0.05, edgecolor='none')
hb = ax.plot(f0[:,2], f0[:,1], marker='o', c='black')

#cb.set_ticks(np.linspace(hb.get_array().min(), hb.get_array().max(), 6))
#cb.set_ticklabels(np.linspace(0, 1., 6))
#cb.set_label('Normalized Count', fontsize=16)
#cb.ax.tick_params(labelsize=16)

#ax.set_title('H2SO4 particle (Number Concentration: '+r' 10$^2 g/m^3$)', fontsize=16)
#ax.set_xlabel('Particle radius ('+'$\mu$'+'m)', fontsize=16)
ax.set_ylabel('Color ratio (520 nm/1022 nm)', fontsize=14)

ax.set_xscale('log')

ax.set_xlabel('Extinction coefficient '+r'$(km^{-1})$'+' @ 1022 nm', fontsize=14)
#ax42.set_ylabel('Normalized Frequency', fontsize=10)
#ax.set_title('Two Hibit model (THM)', fontsize=16)
ax.set_xlim(0.00001,0.1)
ax.set_ylim(0,6.)
#ax.set_xticklabels(fontsize=14)
ax.tick_params(labelsize=16)
#ax.set_xticks([0,0.5,1,1.5])

ax.text(f0[0,2]*0.2,f0[0,1]-0.17, '0.08 '+'$\mu$'+'m\n'+str(f0[0,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[1,2]*0.2,f0[1,1]-0.17, '0.1 '+'$\mu$'+'m\n'+str(f0[1,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[2,2]*0.2,f0[2,1]-0.1, '0.15 '+'$\mu$'+'m\n'+str(f0[2,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[3,2]*0.2,f0[3,1]-0.15, '0.2 '+'$\mu$'+'m\n'+str(f0[3,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[4,2]*0.29,f0[4,1]-0.5, '0.3 '+'$\mu$'+'m\n'+str(f0[4,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[5,2]*0.57,f0[5,1]+0.3, '0.4 '+'$\mu$'+'m\n'+str(f0[5,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[6,2]*0.6,f0[6,1]-0.6, '0.5 '+'$\mu$'+'m\n'+str(f0[6,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[7,2]*0.65,f0[7,1]+0.3, '0.75 '+'$\mu$'+'m\n'+str(f0[7,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[8,2]*0.68,f0[8,1]-0.7, '1.0 '+'$\mu$'+'m\n'+str(f0[8,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[9,2]*0.75,f0[9,1]+0.3, '1.5 '+'$\mu$'+'m\n'+str(f0[9,-1])+r' $cm^{-3}$',color='black',fontsize=9)
ax.text(f0[10,2]*0.5,f0[10,1]-0.6, '3.0 '+'$\mu$'+'m\n'+str(f0[10,-1])+r' $cm^{-3}$',color='black',fontsize=9)

#ax.arrow(f0[4,2]+1,f0[4,1]+0.05,0.51,-0.5, head_width=0.15, head_length=0.15, fc='k', ec='k')
#ax.arrow(f0[5,2]+0.05,f0[7,1]-0.5,0.08,0.25, head_width=0.04, head_length=0.04, fc='k', ec='k')

#ax.text(0.0000005, 0.5, 'density = 0.1 '+r'$cm^{-3}$', fontsize = 15)

#x = np.arange(6)
#y = np.arange(6)

#plt.plot(x, y, c='black', linewidth=0.5)
#ax.text(0.2, 4.5, '(b)', fontsize = 18)

pngname = "./"+"plt_h2so4_ExtnMix"+".pdf"
print("save ", pngname)
plt.savefig(pngname, dpi=100, facecolor='w', edgecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
