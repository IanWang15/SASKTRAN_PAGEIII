import numpy as np
import matplotlib.pyplot as plt

f0 = np.loadtxt("./f_h2so4.txt")

fig = plt.figure(figsize=(6.1,5))
ax=fig.add_subplot(111)
#plt.scatter(filedata[:,1], filedata[:,2], c='navy', alpha=0.05, edgecolor='none')
hb = ax.plot(f0[:,0], f0[:,1], c='black')

#cb.set_ticks(np.linspace(hb.get_array().min(), hb.get_array().max(), 6))
#cb.set_ticklabels(np.linspace(0, 1., 6))
#cb.set_label('Normalized Count', fontsize=16)
#cb.ax.tick_params(labelsize=16)

#0.1 4.027328701422502
#0.2 2.0171640542799922
#0.3 1.3001908652969205
#0.5 0.8744683077968979
#1.0 0.8873748271732115
#1.5 0.936821347373896

loclist = [4.027328701422502, 2.0171640542799922, 1.3001908652969205,0.8744683077968979]
ax.plot([0.5,0.5],[0.,loclist[3]], '--', c='black', linewidth=0.7)
ax.plot([0.1,0.1],[0.,loclist[0]], '--', c='black', linewidth=0.7)
ax.plot([0.2,0.2],[0.,loclist[1]], '--', c='black', linewidth=0.7)
ax.plot([0.3,0.3],[0.,loclist[2]], '--', c='black', linewidth=0.7)
ax.plot([0.5,0.], [loclist[3],loclist[3]],'--', c='black', linewidth=0.7)
ax.plot([0.1,0.],[loclist[0],loclist[0]],'--', c='black', linewidth=0.7)
ax.plot([0.2,0.],[loclist[1],loclist[1]], '--', c='black', linewidth=0.7)
ax.plot([0.3,0.],[loclist[2],loclist[2]], '--', c='black', linewidth=0.7)

ax.text(0.14,loclist[0], '(0.1'+'$\mu$'+'m, 4.03)',color='black',fontsize=12)
ax.text(0.24,loclist[1], '(0.2'+'$\mu$'+'m, 2.02)',color='black',fontsize=12)
ax.text(0.34,loclist[2], '(0.3'+'$\mu$'+'m, 1.30)',color='black',fontsize=12)
ax.text(0.51,loclist[3]+0.1, '(0.5'+'$\mu$'+'m, 0.89)',color='black',fontsize=12)

#ax.set_title('H2SO4 particle (Number Concentration: '+r' 10$^2 g/m^3$)', fontsize=16)
ax.set_xlabel('Particle radius ('+'$\mu$'+'m)', fontsize=14)
ax.set_ylabel('Color ratio (520 nm/1022 nm)', fontsize=14)
#ax42.set_ylabel('Normalized Frequency', fontsize=10)
#ax.set_title('Two Hibit model (THM)', fontsize=16)
ax.set_xlim(0,1.5)
ax.set_ylim(0,5.)
#ax.set_xticklabels(fontsize=14)
ax.tick_params(labelsize=16)
ax.set_xticks([0,0.5,1,1.5])

#x = np.arange(6)
#y = np.arange(6)

#plt.plot(x, y, c='black', linewidth=0.5)
#ax.text(0.2, 4.5, '(b)', fontsize = 18)

pngname = "./"+"plt_h2so4_radis2"+".pdf"
print("save ", pngname)
plt.savefig(pngname, dpi=100, facecolor='w', edgecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
