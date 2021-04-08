import numpy as np
import matplotlib.pyplot as plt

f0 = np.loadtxt("./extiCoff.txt")

fig = plt.figure(figsize=(6.1,5))
ax=fig.add_subplot(111)
#plt.scatter(filedata[:,1], filedata[:,2], c='navy', alpha=0.05, edgecolor='none')
hb = ax.plot(f0[:,0], f0[:,-1], c='black')

#cb.set_ticks(np.linspace(hb.get_array().min(), hb.get_array().max(), 6))
#cb.set_ticklabels(np.linspace(0, 1., 6))
#cb.set_label('Normalized Count', fontsize=16)
#cb.ax.tick_params(labelsize=16)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_title('H2SO4 particle, (0.5 '+'$\mu$'+'m)', fontsize=16)
ax.set_ylabel('Log 10 Extinction Coefficient '+r'$(km^{-1})$', fontsize=16)
ax.set_xlabel('Number Concentration '+r'$(g/m^3)$', fontsize=16)
#ax42.set_ylabel('Normalized Frequency', fontsize=10)
#ax.set_title('Two Hibit model (THM)', fontsize=16)
#ax.set_xlim(0,2.5)
#ax.set_ylim(0,5.)
#ax.set_yticklabels(['',0,'',0.05,'',0.1,'',0.15])
ax.tick_params(labelsize=16)

#x = np.arange(6)
#y = np.arange(6)

#plt.plot(x, y, c='black', linewidth=0.5)
#ax.text(0.2, 4.5, '(b)', fontsize = 18)

pngname = "./"+"plt_h2so4_extn"+".pdf"
print("save ", pngname)
plt.savefig(pngname, dpi=100, facecolor='w', edgecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
