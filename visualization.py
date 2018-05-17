import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text
from matplotlib.patches import FancyArrow, Rectangle
import matplotlib
import argparse
import os
ip=1



parser = argparse.ArgumentParser(description='PiPred')

parser.add_argument('-i',
                    help='Path to output file from PiPred',
                    required=True,
					metavar='DIR',
					default = '.')

parser.add_argument('-m',
                    help='Plot mode: heatmap, bar, curve',
					default = 'curve')
args = parser.parse_args()

def read_output(out_file):
	

	data = pd.read_csv(out_file,sep ='\s+',index_col=0)
	return data,list(data.index.values),out_file.split('.')[0]


def gen_heatmap(data,savefig=False):
	heatmap = sns.heatmap(data)
	plt.show()
	


def gen_probability_graph(prob,seq,savefig=False):
	
	def pi_plot(ind,value):
		ax=fig.add_subplot(1,1,1)
	
		ax.set_ylim([0,1])

	
		ax.tick_params(axis='x',          
		which='major',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',         # ticks along the top edge are off
		labelbottom='on')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		for cv in range(ind):
			if np.argmax(prob[cv])==2:
				pi = ax.bar(cv,I[cv+value*ind],w,color=(0,0,1),picker=True)
			else:
				pi = ax.bar(cv,I[cv+value*ind],w,color=(0.6,0.6,1),picker=True)
		

		ax.set_xticklabels(seq[value*ind:(value+1)*ind])
		ax.set_xticks(range(value*ind,(value+1)*ind))

		for label in ax.get_xticklabels():  # make the xtick labels pickable
			label.set_picker(True)
		art = ax.arrow(ind-2,0.96,1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='next')
		art1 = art = ax.arrow(ind-2,0.96,-1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='back')	

		fig.text(0.08, 0.5, 'probability', ha='center', va='center', rotation='vertical')
		fig.text(0.5, 0.9, 'residues {} - {}'.format(value*ind+1,(value+1)*ind), ha='center', va='center', rotation='horizontal')
	
		ax1 = fig.add_axes((0,0,1,0.06))
		ax1.spines['top'].set_visible(False)
		ax.set_picker(True)
	nb = 0


	w=0.8
	
	seq_len = len(seq)
	ind=100
	tot = seq_len//ind+1
	rest = tot*ind
	print(prob.shape,len(seq))
	prob = np.concatenate((prob,np.zeros(((tot)*ind-seq_len,4))))
	
	for i in range((tot)*ind-seq_len):
		seq.append(' ')
	prt =seq[:]
	E,H,I,C=prob[:,0],prob[:,1],prob[:,2],prob[:,3]
	fig = plt.figure(figsize=(30,15))
	pi_plot(ind,0)
	
	
				
	def onpick(event):
	
		global ip
		patch =event.artist
		if patch:
			pop=0

			ax1 = fig.add_axes((0,0,1,0.06))
			ax1.spines['top'].set_visible(False)
	
		
			if isinstance(patch,Text):
				ax1.clear()
				ax1.spines['top'].set_visible(False)
				wt = patch.get_position()[0]
				wr= patch.get_text()
				pb,ph,pi,pc=round(prob[(ip-1)*ind+wt][0],4),round(prob[(ip-1)*ind+wt][1],4),round(prob[(ip-1)*ind+wt][2],4),round(prob[(ip-1)*ind+wt][3],4)
				ax1.patch.set_alpha(0.1)
				ax1.patch.set_color('b')
		
				
				nap = 'residue nr {},aa_name {}, beta: {} alpha: {} pi:{} coil: {}'.format((ip-1)*ind+wt+1,wr,str(pb),str(ph),str(pi),str(pc))
				ax1.text(0.5, 0.2, nap, ha='center', va='center')

			elif isinstance(patch,FancyArrow) and patch._label == 'next':
				if ip<tot:
					fig.clear()
					ax1 = fig.add_axes((0,0,1,0.06))
					ax1.spines['top'].set_visible(False)
					ax = fig.add_subplot(1,1,1)
					ax.set_picker(True)
					ax.set_ylim([0,1])
			
				
					ax.tick_params(axis='x',          
					which='major',      # both major and minor ticks are affected
					bottom='off',      # ticks along the bottom edge are off
					top='off',         # ticks along the top edge are off
					labelbottom='on')
					ax.spines['right'].set_visible(False)
					ax.spines['top'].set_visible(False)
					for cv in range(ind):
						if np.argmax(prob[ip*ind+cv])==2:
							pi = ax.bar(cv,I[ip*ind+cv],w,color=(0,0,1),picker=True)
						else:
							pi = ax.bar(cv,I[ip*ind+cv],w,color=(0.6,0.6,1),picker=True)
				
			
				
					ax.set_xticklabels(seq[ip*ind:ip*ind+ind])
					ax.set_xticks(range(ind))

					for label in ax.get_xticklabels():  # make the xtick labels pickable
						label.set_picker(True)
					#ax.clear()
					art = ax.arrow(ind-2,0.96,1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='next')
					art1 = art = ax.arrow(ind-2,0.96,-1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='back')
					fig.text(0.08, 0.5, 'probability', ha='center', va='center', rotation='vertical')
					if ip+1 == tot:
						fig.text(0.5, 0.9, 'residues {} - {}'.format(ip*ind+1,seq_len), ha='center', va='center', rotation='horizontal')
					else:
						fig.text(0.5, 0.9, 'residues {} - {}'.format(ip*ind+1,ip*ind+ind), ha='center', va='center', rotation='horizontal')
					ip+=1

			elif isinstance(patch, FancyArrow) and patch._label == 'back':
				if ip>1:
					ip-=1
					fig.clear()
					ax1 = fig.add_axes((0,0,1,0.06))
					ax1.spines['top'].set_visible(False)
					ax = fig.add_subplot(1,1,1)
					ax.set_picker(True)
					ax.set_ylim([0,1])

				
					ax.tick_params(axis='x',          
					which='major',     
					bottom='off',      
					top='off',         
					labelbottom='on')
					ax.spines['right'].set_visible(False)
					ax.spines['top'].set_visible(False)
				
					for cv in range(ind):
						if np.argmax(prob[(ip-1)*ind+cv])==2:
							pi = ax.bar(cv,I[(ip-1)*ind+cv],w,color=(0,0,1),picker=True)
						else:
							pi = ax.bar(cv,I[(ip-1)*ind+cv],w,color=(0.6,0.6,1),picker=True)
				
			
				
					ax.set_xticklabels(seq[(ip-1)*ind:(ip)*ind])
					ax.set_xticks(range(ind))
					for label in ax.get_xticklabels():  
						label.set_picker(True)
					art = ax.arrow(ind-2,0.96,1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='next')
					art1 = art = ax.arrow(ind-2,0.96,-1.5,0,width=0.016,head_width=0.034,head_length=0.24,picker=True,label='back')
					
					fig.text(0.08, 0.5, 'probability', ha='center', va='center', rotation='vertical')
					fig.text(0.5, 0.9, 'residues {} - {}'.format((ip-1)*ind+1,(ip)*ind), ha='center', va='center', rotation='horizontal')
					#ax.add_patch(rect)
				
				
			else:
			
				ax1.clear()
			
			fig.canvas.draw_idle()
		
		
	fig.canvas.mpl_connect('pick_event',onpick)
	plt.show()

def plotCurve(prob,title,savefig=False,wyb = ['I'],figsize = (10,5)):
		

		fig = plt.figure(figsize=figsize)
		dl = len(prob)
		lab = {'E':0,'H':1,'I':2,'C':3}
		ax = fig.add_subplot(1,1,1)
		ax.set_ylim([0,1])
		ax.set_title(title)
		ax.set_xlabel('residue_number')
		ax.set_ylabel('pi_probability')
		for i in wyb:
			ax.plot(range(dl),prob[:,lab[i]],label=i)
		ax.legend()
		if savefig:
			plt.savefig('prob_dist1',dpi=300)
		plt.show()

	


if __name__ == '__main__':
	file_name = args.i
	mode = args.m
	prob,seq,title = read_output(file_name)
	#prob = np.array(prob)
	if mode == 'curve':
		plotCurve(np.array(prob),title,savefig = True)
	elif mode == 'heatmap':
		gen_heatmap(prob)
	elif mode =='bar':
		gen_probability_graph(prob,seq)




