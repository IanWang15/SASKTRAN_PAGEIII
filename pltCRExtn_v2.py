import numpy as np
import matplotlib.pyplot as plt

f0 = np.loadtxt("./f_h2so4.txt")

fig = plt.figure(figsize=(6.1,5))
ax=fig.add_subplot(111)
#plt.scatter(filedata[:,1], filedata[:,2], c='navy', alpha=0.05, edgecolor='none')
hb = ax.plot(f0[:,-1], f0[:,1], c='black')

#cb.set_ticks(np.linspace(hb.get_array().min(), hb.get_array().max(), 6))
#cb.set_ticklabels(np.linspace(0, 1., 6))
#cb.set_label('Normalized Count', fontsize=16)
#cb.ax.tick_params(labelsize=16)

loclist = [4.027328701422502, 2.0171640542799922, 1.3001908652969205,0.8744683077968979]
ax.plot([f0[1,-1],f0[1,-1]],[0.,f0[1,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[3,-1],f0[3,-1]],[0.,f0[3,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[5,-1],f0[5,-1]],[0.,f0[5,1]], '--', c='black', linewidth=0.7)
#ax.plot([f0[7,-1],f0[7,-1]],[0.,f0[7,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[9,-1],f0[9,-1]],[0.,f0[9,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[19,-1],f0[19,-1]],[0.,f0[19,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[29,-1],f0[29,-1]],[0.,f0[29,1]], '--', c='black', linewidth=0.7)

ax.plot([f0[1,-1],0],[f0[1,1],f0[1,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[3,-1],0],[f0[3,1],f0[3,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[5,-1],0],[f0[5,1],f0[5,1]], '--', c='black', linewidth=0.7)
#ax.plot([f0[7,-1],0],[f0[7,1],f0[7,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[9,-1],0],[f0[9,1],f0[9,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[19,-1],0],[f0[19,1],f0[19,1]], '--', c='black', linewidth=0.7)
ax.plot([f0[29,-1],0],[f0[29,1],f0[29,1]], '--', c='black', linewidth=0.7)

ax.text(f0[1,-1]+0.0001,f0[1,1], 'r = 0.1 '+'$\mu$'+'m',color='black',fontsize=12)
ax.text(f0[3,-1]+0.002,f0[3,1]+0.02, 'r = 0.2 '+'$\mu$'+'m',color='black',fontsize=12)
ax.text(f0[5,-1],f0[5,1]+0.05, 'r = 0.3 '+'$\mu$'+'m',color='black',fontsize=12)
ax.text(f0[4,-1]-0.0019,f0[7,1]-0.7, 'r = 0.5 '+'$\mu$'+'m',color='black',fontsize=12)
ax.text(f0[7,-1]+0.01,f0[7,1]+0.05, 'r = 1.0 '+'$\mu$'+'m',color='black',fontsize=12)
ax.text(f0[4,-1]+0.2,f0[4,1]+0.2, 'r = 1.5 '+'$\mu$'+'m',color='black',fontsize=12)
ax.arrow(f0[4,-1]+1,f0[4,1]+0.05,0.51,-0.5, head_width=0.15, head_length=0.15, fc='k', ec='k')
ax.arrow(f0[5,-1]+0.05,f0[7,1]-0.5,0.08,0.25, head_width=0.04, head_length=0.04, fc='k', ec='k')

#ax.set_title('H2SO4 particle (Number Concentration: '+r' 10$^2 g/m^3$)', fontsize=16)
#ax.set_xlabel('Particle radius ('+'$\mu$'+'m)', fontsize=16)
ax.set_ylabel('Color ratio (520 nm/1022 nm)', fontsize=14)

ax.set_xscale('log')

ax.set_xlabel('Extinction coefficient '+r'$(km^{-1})$'+' @ 1022 nm', fontsize=14)
#ax42.set_ylabel('Normalized Frequency', fontsize=10)
#ax.set_title('Two Hibit model (THM)', fontsize=16)
ax.set_xlim(0.0001,2.)
ax.set_ylim(0,5.)
#ax.set_xticklabels(fontsize=14)
ax.tick_params(labelsize=16)
#ax.set_xticks([0,0.5,1,1.5])

#x = np.arange(6)
#y = np.arange(6)

#plt.plot(x, y, c='black', linewidth=0.5)
#ax.text(0.2, 4.5, '(b)', fontsize = 18)

pngname = "./"+"plt_h2so4_Extn"+".pdf"
print("save ", pngname)
plt.savefig(pngname, dpi=100, facecolor='w', edgecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
