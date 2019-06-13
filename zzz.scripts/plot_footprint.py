import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

fileName = sys.argv[1]
cutoff = float(sys.argv[2])
maxres = int(sys.argv[3])

ZINClist = open(fileName,'r')
ZINCids = ZINClist.readlines()
ZINClist.close()

counter = []

temp = open(ZINCids[0].strip()+'.txt','r')
templines = temp.readlines()
temp.close()

for templine in templines:
	templinesplit = templine.split()
	if (len(templinesplit) == 8):
		if (templinesplit[0] != 'resname'):
			counter.append(0)

for ZINCid in ZINCids:
	footprint = open(ZINCid.strip()+'.txt','r')
	lines = footprint.readlines()
	footprint.close()

	vdw_ref = []
	vdw_pose = []
	es_ref = []
	es_pose = []
	resname = []
	resid = []

	for line in lines:
		linesplit = line.split()
		if (len(linesplit) == 8):
			if (linesplit[0] != 'resname'):
				resname.append(linesplit[0])
				resid.append(int(linesplit[1]))
				vdw_ref.append(float(linesplit[2]))
				es_ref.append(float(linesplit[3]))
				vdw_pose.append(float(linesplit[5]))
				es_pose.append(float(linesplit[6]))

	for i in range(len(resid)):
		if (vdw_ref[i] < -cutoff or vdw_pose[i] < -cutoff or es_ref[i] < -cutoff or es_pose[i] < -cutoff or vdw_ref[i] > cutoff or vdw_pose[i] > cutoff or es_ref[i] > cutoff or es_pose[i] > cutoff):
			counter[i] += 1


b = counter[:]
indice = []
for i in range(maxres):
	tmp = max(b)
	temp = counter.index(tmp)
	while temp+1 in indice:
		temp = counter[temp+1:len(counter)].index(tmp) + temp + 1
	indice.append(temp+1)
	b.remove(tmp)
indice.sort()

resname_selected = []
j = 0

for i in range(len(resid)):
	if (resid[i] == indice[j]):
		resname_selected.append(resname[i])
		if j < len(indice)-1:
			j += 1

resname_selected.append('REMAIN')

for ZINCid in ZINCids:
	footprint = open(ZINCid.strip()+'.txt','r')
        lines = footprint.readlines()
        footprint.close()

	vdw_ref = []
        vdw_pose = []
        es_ref = []
        es_pose = []
        resname = []
        resid = []

        for line in lines:
                linesplit = line.split()
		if (len(linesplit) == 3):
			if (linesplit[1] == 'vdw_fp:'):
				vdw_score = 'd = '+linesplit[2]
			if (linesplit[1] ==  'es_fp:'):
				es_score = 'd = '+linesplit[2]
                if (len(linesplit) == 8):
                        if (linesplit[0] != 'resname'):
                                resname.append(linesplit[0])
                                resid.append(int(linesplit[1]))
                                vdw_ref.append(float(linesplit[2]))
                                es_ref.append(float(linesplit[3]))
                                vdw_pose.append(float(linesplit[5]))
                                es_pose.append(float(linesplit[6]))

	resid_selected = []
	vdw_ref_selected = []
	vdw_pose_selected = []
	es_ref_selected = []
	es_pose_selected = []
	vdw_ref_left = 0
	vdw_pose_left = 0
	es_ref_left = 0
	es_pose_left = 0
	j = 0

	for i in range(len(resname)):
		if (resname[i] == resname_selected[j]):
			resid_selected.append(resid[i])
			vdw_ref_selected.append(vdw_ref[i])
			vdw_pose_selected.append(vdw_pose[i])
			es_ref_selected.append(es_ref[i])
			es_pose_selected.append(es_pose[i])
			if j < len(resname_selected)-1:
				j += 1
		else:
			vdw_ref_left = vdw_ref_left + vdw_ref[i]
			vdw_pose_left = vdw_pose_left + vdw_pose[i]
			es_ref_left = es_ref_left + es_ref[i]
	                es_pose_left = es_pose_left + es_pose[i]


	vdw_ref_selected.append(vdw_ref_left)
	vdw_pose_selected.append(vdw_pose_left)
	es_ref_selected.append(es_ref_left)
	es_pose_selected.append(es_pose_left)
	
	residue = []
	
	for i in range(len(resname_selected)):
		residue.append(i)
	
	# print (sum(vdw_pose_selected) - sum(vdw_pose)) + (sum(vdw_ref_selected) - sum(vdw_ref)) + (sum(es_pose_selected) - sum(es_pose)) + (sum(es_ref_selected) - sum(es_ref))
		
	fig = plt.figure(figsize=(12, 11))
	ax1 = fig.add_subplot(2,1,1)
	ax1.set_title(ZINCid.strip())
	plt.plot(residue, vdw_ref_selected, 'b', linewidth=3)
	plt.plot(residue, vdw_pose_selected, 'r', linewidth=3)
	ax1.set_ylabel('VDW Energy')
	ax1.set_ylim(-10, 5)
	ax1.set_xlim(0, len(resname_selected))
	ax1.xaxis.set_major_locator(MultipleLocator(1))
	ax1.xaxis.set_major_formatter(FormatStrFormatter('%s'))
	ax1.set_xticks(residue)
	ax1.xaxis.grid(which='major', color='black', linestyle='solid')
	ax1.set_xticklabels(resname_selected, rotation=90)
	ax1.legend(['Reference', 'Pose'])
	ax1.annotate(vdw_score, xy=(40,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})

	ax2 = fig.add_subplot(2,1,2)
	plt.plot(residue, es_ref_selected, 'b', linewidth=3)
	plt.plot(residue, es_pose_selected, 'r', linewidth=3)
	ax2.set_ylabel('ES Energy')
	ax2.set_ylim(-10, 5)
	ax2.set_xlim(0, len(resname_selected))
	ax2.xaxis.set_major_locator(MultipleLocator(1))
	ax2.xaxis.set_major_formatter(FormatStrFormatter('%s'))
	ax2.set_xticks(residue)
	ax2.xaxis.grid(which='major', color='black', linestyle='solid')
	ax2.set_xticklabels(resname_selected, rotation=90)
	ax2.legend()
	ax2.annotate(es_score, xy=(40,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
	
	plt.savefig(ZINCid.strip()+'.pdf')

