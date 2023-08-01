# This is a python script to check the rigidity of a 2D frictionless 
# disk packing. I will write this in python as if this is a pseudo code
# that works.

# HERE I DEFINE SOME VARIABLES AND WHAT THEY DO !
#
# pebbles is a N X 2 array : stores 2 pebbles for each site : pebbles[n][0] would store the ID
# of the site with which n shares an edge, and that edge is covered using a pebble from n. if 
# pebbles[n][0]=None then that pebble is free ! similarly for pebbles[n][1] which stores info 
# about the other pebble
#
# path is an N array : path(n) will store which site one should proceed from site n. path(n)=-1 
# means the path is dead end. ( No directed graph from n)
#
# seen is an N array : seen(v)= True/False : stores if we have been on the site v ( looking for free 
# pebble ?? !)




#####################NEW UPDATE##################################
##
##			THIS CODE AS OF NOW IMPLEMENTS GENERAL (K,L) game
##					(untested )
##
#################################################################
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import collections  as mc
import math
import random
sys.setrecursionlimit(10**6)
def Find_Pebble(v,seen,path,Pebbles):
    # We have seen the pebble 
    seen[v]=True
    #print(K)  
    # We have the path reset : v does not point to anyone
    path[v]=-1
    # Check if any of the K pebbles are unassigned/free
    for i in range(0,K):
        if(Pebbles[v][i]==None):
            return True # Return that a free pebble is found 
    #if(Pebbles[v][1]==None):
    #    return True # the other pebble is free and hence return
    for i in range(0,K):
        X=Pebbles[v][i] # X is the other node with which v has an edge and that edge is covered with the i'th pebble of v
        if( seen[X] == False ):
            # so the directed path from v  points towards X
            path[v]=X
            if(Find_Pebble(X,seen,path,Pebbles)): # go look for a free pebble from X
                return True 
 ## Y=Pebbles[v][1] 
 ## if( seen[Y] == False ): 
 ##     path[v]=Y
 ##     if(Find_Pebble(Y,seen,path,Pebbles)): # go look for a free pebble from Y
 ##         return True 
    return False

def Rearrange_Pebbles(v,v_,path,Pebbles):
    # This function should rearrange pebbles after a pebble has been found 
    # If one follows path[] one will reach a node with a free pebble !!
    # v_ is the associated vertice with the edge we are considering. 
    for i in range (0,K):
        if(Pebbles[v][i]==None): # Let us just check if there is a free pebble and move on !
           Pebbles[v][i]=v_
           return 
    #if(Pebbles[v][1]==None): # Just do it for the other pebble 
    #    Pebbles[v][1]=v_
    #    return 
    v_copy=v # We will have to store this guy because, after traversing the path and changing
             # all the directions (pebble arrangements), we need to change for the first guy 
    while(path[v]!=-1): # Start following the path 
        X=path[v] # this is where it point to !
        if(path[X] != -1 ): # If the node to which v points to is not the end of path 
            W=path[X]
            for i in range(0,K):
                if(Pebbles[X][i]==W): # find the pebble which is used to cover w ( from x)
                    Pebbles[X][i]=v
                    break
            #if(Pebbles[X][1]==W): # find the pebble which is used to cover w ( from x)
            #    Pebbles[X][1]=v
        else: # Thus the path[X] is the last in the path , this guy should have a free pebble
            for i in range (0,K):
                if(Pebbles[X][i]==None): # this is free pebble
                    Pebbles[X][i]=v
                    break
            #elif(Pebbles[X][1]==None): # check other
            #    Pebbles[X][1]=v
        v=X
    # Now that we have reversed the directions of the directed path we can change the pebble 
    # arrangemment at the very first node of the directed path 
    for i in range(0,K):
        if(Pebbles[v_copy][i]==path[v_copy]):
            Pebbles[v_copy][i]=v_
            return 
    #if(Pebbles[v_copy][1]==path[v_copy]):
    #    Pebbles[v_copy][1]=v_
    #    return 

def PebbleGame(e,Pebbles,seen,path,independent_edges,redundant_edges):
    #for e in G:
    Pebbles_copy=Pebbles.copy()
    flag=1 # lets assume that the new edge is independent 
    # we will work with a copy of the current pebble arrangement
    for i in range (0,L+1):
        #print i 
        seen=np.full(N,False)
        path=np.full(N,-1)
        #found = 
        #print path,seen,found
        if(Find_Pebble(e[0],seen,path,Pebbles_copy)):
            Rearrange_Pebbles(e[0],e[1],path,Pebbles_copy)
            seen=np.full(N,False)
            path=np.full(N,-1)
        elif(seen[e[1]]==False and Find_Pebble(e[1],seen,path,Pebbles_copy) ):
            Rearrange_Pebbles(e[1],e[0],path,Pebbles_copy)
            seen=np.full(N,False)
            path=np.full(N,-1)
        else:
            if(i == L): 
                #print "redundant bond"
                sum=0
                for i in range (0,N):
                    if(seen[i]):
                        #sum=sum+1
                        over_const_sites[i]=True
                over_const_per_red.append(over_const_sites)
                redundant_edges.append(e)
                #print sum
         ##     over_constrained.append(e[0])
         ##     print "begin"
         ##     print e[0]
         ##     next=path[e[0]]
         ##     while( next != -1 ):
         ##         print next
         ##         #over_constrained.append([this,next])
         ##         over_constrained.append(next)
         ##         #this=next
         ##         next=path[next]
         ##     print "end"
         ##     if(seen[e[1]] == False):
         ##         print "anytime here?"
         ##         #this=e[1]
         ##         over_constrained.append(e[1])
         ##         next=path[e[1]]
         ##         while( next != -1 ):
         ##                 #over_constrained.append([this,next])
         ##                 over_constrained.append(next)
         ##                 #this=next
         ##                 next=path[next]
            if(i < L): 
                print("should not happen")
            flag=0 # well the edge is not independent and hence redundant
            break 
        #print Pebbles_copy
    # Now we need to remove the three extra edges and keep one
    # We will do that here, using the orginal orientation of pebbles, before quadruple busineess 
    if(flag):
        #print path,seen,found
        independent_edges.append(e)
        if(Find_Pebble(e[0],seen,path,Pebbles)):
            Rearrange_Pebbles(e[0],e[1],path,Pebbles)
        elif(seen[e[1]]==False and Find_Pebble(e[1],seen,path,Pebbles)):
            Rearrange_Pebbles(e[1],e[0],path,Pebbles)
        else :
            print("error")
        #print Pebbles
    free_pebbles=0
    for i,j in enumerate(Pebbles):
        for k in range(0,K):
            if(j[k]==None):
                free_pebbles=free_pebbles+1
        #if(j[1]==None):
        #    free_pebbles=free_pebbles+1
    #print free_pebbles
    return free_pebbles
    #print e,len(independent_edges), len(redundant_edges),free_pebbles



## MAIN PART OF THE CODE STARTS HERE 

K=3
L=2

N=int(sys.argv[2])
Neibhours= [ [] for i in range(N) ]
Pebbles=np.full((N,K),None)
seen=np.full(N,False)
over_constrained=[]
over_const_per_red=[]
over_const_sites=np.full(N,False)
path=np.full(N,-1)
independent_edges=[]
redundant_edges=[]
#over_constrained=[]
free_pebbles=0
rigid=[] # these are empty
floppy=[]
marked=np.full(N,0) # NONE OF THE SITES HAS BEEN SEEN
EDGES=np.zeros(shape=(N,N))
I_EDGES=np.zeros(shape=(N,N))
RIGID_EDGES=np.full((N,N),-1)
OVER_CONST_EDGES=np.full((N,N),-1)
#for e in edge_set:
#    PebbleGame(e,Pebbles,seen,path,independent_edges,redundant_edges)
F=open(sys.argv[1],"r")
edges=F.readlines()
avg_z = len(edges)*1./N
for f in range (0,len(edges)):
    e=[]
    e.append(int(edges[f].split('\t')[0]))
    e.append(int((edges[f].split('\t')[1]).split('\n')[0]))
    free_pebbles=PebbleGame(e,Pebbles,seen,path,independent_edges,redundant_edges)
    Neibhours[e[0]].append(e[1])
    Neibhours[e[1]].append(e[0])
    EDGES[e[0]][e[1]]=1
    EDGES[e[1]][e[0]]=1
print(free_pebbles)
#exit()
# At this point all edges are identified as either independent or redundant. 
# Now these independent edges need to be decomposed into rigid cluster. 
# The algorithm is given 
#   Consider an independent bond, which is not labeled ( as a rigid or floppy site ??!! )
#
#   #break
#print Neibhours
for os1 in range (0,N):
    for os2 in range (0,N):
        if(over_const_sites[os1] and over_const_sites[os2]):
            if(EDGES[os1][os2]):
                OVER_CONST_EDGES[os1][os2]=1
                OVER_CONST_EDGES[os2][os1]=1
                over_constrained.append([os1,os2])
            
        
for label,ie in enumerate(independent_edges):
    I_EDGES[ie[0]][ie[1]]=1
    I_EDGES[ie[1]][ie[0]]=1
rigid_clusters = [ [] for _ in independent_edges ]
max_size=-1
max_size_C=-1
for label,ie in enumerate(independent_edges):
    #print label,ie
    rigid=[] # these are empty
    floppy=[]
    marked=np.full(N,0) # NONE OF THE SITES HAS BEEN SEEN
    Pebbles_copy=Pebbles.copy()
    if(RIGID_EDGES[ie[0]][ie[1]]== -1):
        for i in range (0,L):
            #print i 
            seen=np.full(N,False)
            path=np.full(N,-1)
            #found = 
            #print path,seen,found
            if(Find_Pebble(ie[0],seen,path,Pebbles_copy)):
                Rearrange_Pebbles(ie[0],ie[1],path,Pebbles_copy)
                seen=np.full(N,False)
                path=np.full(N,-1)
            elif(seen[e[1]]==False and Find_Pebble(ie[1],seen,path,Pebbles_copy) ):
                Rearrange_Pebbles(ie[1],ie[0],path,Pebbles_copy)
                seen=np.full(N,False)
                path=np.full(N,-1)
            else:
                #if(i == 3): 
                    #print "redundant bond"
                #    redundant_edges.append(e)
                if(i < L): 
                    print("th should not happen")
                flag=0 # well the edge is not independent and hence redundant
                break 
        #Now three pebbles have been found and locked up with these ie[0] and ie[1].
        rigid.append(ie[0])
        rigid.append(ie[1]) # Added these guys to rigid sites.
        marked[ie[0]]=1
        marked[ie[1]]=1
        # Now we will loop through rigid sites 
        for r in rigid:
            for nn in Neibhours[r]:
                #print r,nn,marked[nn]
                if(marked[nn]==0):
                    marked[nn]=1
                    seen=np.full(N,False)
                    path=np.full(N,-1)
                    if(Find_Pebble(nn,seen,path,Pebbles_copy)): # have to use pebbles copy because thats where the pebbles are locked up with ie
                        # Ha we have a free pebble, this site is floppy wrt to ie[0] and ie[1] and so is everything in the path 
                        floppy.append(nn) # add to floppy
                        next=path[nn]
                        while( next != -1 ):
                            if(marked[next]==0):
                                floppy.append(next)
                                marked[next]=1
                            next=path[next]
                    else:
                        # Free pebble is not found and hence the site is rigid so is everything in the path 
                        #print nn,path
                        rigid.append(nn) # add to rigid
                        next=path[nn]
                        while( next != -1 ): # ad everythin in path
                            if(marked[next]==0):
                                #print "here ",next
                                rigid.append(next)
                                marked[next]=1
                            next=path[next]
        for v in rigid:
            for v_ in rigid:
                if(EDGES[v][v_] and v > v_):
                    RIGID_EDGES[v ][v_]=label
                    RIGID_EDGES[v_][v ]=label
                    rigid_clusters[label].append([v,v_])
        if(len(rigid)):
            if(len(rigid)>max_size):
                max_size=len(rigid)
                max_size_C=len(rigid_clusters[label])
               
        #print rigid,label
with open('redundant_edges_'+(sys.argv[3]),"w") as filename:
    for e in redundant_edges:
        filename.write(str(e)+"\n")
with open('independent_edges'+(sys.argv[3]),"w") as filename:
    for e in independent_edges:
        filename.write(str(e)+"\n")
with open(sys.argv[4],"a") as filename:
    filename.write(str(avg_z)+"\t"+sys.argv[3]+"\t"+str(len(independent_edges))+"\t"+str(len(redundant_edges))+"\t"+str(free_pebbles)+"\t"+str(max_size)+"\n")
##########################################################################################################
############################CODE TO DRAW CONFIGS##########################################################
#################################BEGIN####################################################################
##########################################################################################################

rigid_cluster_dist = open('rigid_cluster_size_%s.dat' %(sys.argv[3]),"w")
for r in rigid_clusters:
    rigid_cluster_dist.write(str(len(r))+"\n")
#%matplotlib notebook
#fig1, ax1 = plt.subplots(1, 1)
#ax1.set_aspect(1)
#edges_to_draw=[]
#colors=[]
#filename=open('rigid_config.xyz',"w") 
#E=np.loadtxt('edges')
#sum1=0
##for r in rigid_clusters:
##    sum1 = sum1 + len(r)
##print(sum1)
#c=0
#for r in rigid_clusters:
#    if(len(r)==max_size_C):
###      ##size=len(r)/(2*N-3)
##        c1=random.uniform(0,1)
##        c2=random.uniform(0,1)
##        c3=random.uniform(0,1)
##        color=(c1,c2,c3,1)
#        for e in r:
#            #print(e)
#            for e_ in E:
#                if(e[0]==int(e_[0]) and e[1]==int(e_[1]) or e[1]==int(e_[0]) and e[0]==int(e_[1])):
#                    #print(e,e_[2:6])
#                    edges_to_draw.append([(e_[2],e_[3]),(e_[4],e_[5])])
#                    c=c+1
            #print(e)
##            #edges_to_draw.append([(CONFIG[e[0]][2],CONFIG[e[0]][3]),(CONFIG[e[1]][2],CONFIG[e[1]][3])])
##            filename.write(str(CONFIG[e[0]][0])+"\t"+str(CONFIG[e[0]][1])+"\t"+str(CONFIG[e[0]][2])+"\t"+str(CONFIG[e[0]][3])+"\t"+str(CONFIG[e[0]][4])+"\n")
##            filename.write(str(CONFIG[e[1]][0])+"\t"+str(CONFIG[e[1]][1])+"\t"+str(CONFIG[e[1]][2])+"\t"+str(CONFIG[e[1]][3])+"\t"+str(CONFIG[e[1]][4])+"\n")
##            x1=CONFIG[e[0]][2]
##            y1=CONFIG[e[0]][3]
##            x2=CONFIG[e[1]][2]
##            y2=CONFIG[e[1]][3]
##            delta_rx=x1-x2
##            delta_ry=y1-y2
##
##            delta_rx=(delta_rx-(tilt*round(delta_ry/BOX)));
##            delta_rx=(delta_rx-(BOX*round(delta_rx/BOX)));
##            delta_ry=(delta_ry-(BOX*round(delta_ry/BOX)));
#             #edges_to_draw.append([(x1,y1),(x1-delta_rx,y1-delta_ry)])
##            colors.append(color)
##            #print e
##            #print CONFIG[e[0],2:4],CONFIG[e[1],2:4]
#        colors.append(color)
#print(c)
#lc = mc.LineCollection(edges_to_draw,colors=colors)
#ax1.add_collection(lc)
##plt.plot()
##plt.show()
#plt.savefig('network.pdf')
#CONFIG=np.loadtxt(sys.argv[5])
#
#BOX=2*float(sys.argv[4])
#gamma=sys.argv[3]
#print(BOX,gamma)
#tilt=float(gamma)*float(BOX)
##print r
##fig, ax = plt.subplots(1, 1)
##plt.ioff()
##ax.set_aspect(1)
##for f in CONFIG:
##    f[2]=(f[2]-(tilt*round(f[3]/BOX)));
##    f[2]=(f[2]-(BOX*round(f[2]/BOX)));
##    f[3]=(f[3]-(BOX*round(f[3]/BOX)));
##    if(f[1]==1):
##        circle=plt.Circle((f[2],f[3]),0.5)
##    if(f[1]==2):
##        circle=plt.Circle((f[2],f[3]),0.7)
##    ax.add_artist(circle)
##ax.set_xlim([min(CONFIG[:,2])-10,max(CONFIG[:,2])+10])
##ax.set_ylim([min(CONFIG[:,3])-10,max(CONFIG[:,3])+10])
#edges_to_draw=[]
#colors=[]
#filename=open('rigid_config.xyz',"w") 
rigid_cluster=open('rigid_edges_%s.xyz'%(sys.argv[3]),"w") 
over_c_f=open('over_constrained_%s.xyz'%(sys.argv[3]),"w") 
#rigid_cluster_dist = open('rigid_cluster_size_%s.dat' %(gamma),"w")
for r in rigid_clusters:
#    rigid_cluster_dist.write(str(len(r))+"\n")
    if(len(r)==max_size_C and len(r)>2):
#      ##size=len(r)/(2*N-3)
#        c1=1.#random.uniform(0,1)
#        c2=0.#random.uniform(0,1)
#        c3=0.#random.uniform(0,1)
#        color=(c1,c2,c3,1)
        for e in r:
#            #edges_to_draw.append([(CONFIG[e[0]][2],CONFIG[e[0]][3]),(CONFIG[e[1]][2],CONFIG[e[1]][3])])
            rigid_cluster.write(str(e[0])+"\t"+str(e[1])+"\n")
#            #filename.write(str(CONFIG[e[0]][0])+"\t"+str(CONFIG[e[0]][1])+"\t"+str(CONFIG[e[0]][2])+"\t"+str(CONFIG[e[0]][3])+"\t"+str(CONFIG[e[0]][4])+"\n")
#            #filename.write(str(CONFIG[e[1]][0])+"\t"+str(CONFIG[e[1]][1])+"\t"+str(CONFIG[e[1]][2])+"\t"+str(CONFIG[e[1]][3])+"\t"+str(CONFIG[e[1]][4])+"\n")
##            x1=CONFIG[e[0]][2]
##            y1=CONFIG[e[0]][3]
##            x2=CONFIG[e[1]][2]
##            y2=CONFIG[e[1]][3]
##            delta_rx=x1-x2
##            delta_ry=y1-y2
##
##            delta_rx=(delta_rx-(tilt*round(delta_ry/BOX)));
##            delta_rx=(delta_rx-(BOX*round(delta_rx/BOX)));
##            delta_ry=(delta_ry-(BOX*round(delta_ry/BOX)));
##            edges_to_draw.append([(x1,y1),(x1-delta_rx,y1-delta_ry)])
##            colors.append(color)
##            #print e
##            #print CONFIG[e[0],2:4],CONFIG[e[1],2:4]
#color=(0.,0.,1.,1)
for oe in over_constrained:
     over_c_f.write(str(oe[0])+"\t"+str(oe[1])+"\n")
     #print(oe)
     #print re
#     x1=CONFIG[oe[0]][2]
#     y1=CONFIG[oe[0]][3]
#     x2=CONFIG[oe[1]][2]
#     y2=CONFIG[oe[1]][3]
#     delta_rx=x1-x2
#     delta_ry=y1-y2
# 
#     delta_rx=(delta_rx-(tilt*round(delta_ry/BOX)));
#     delta_rx=(delta_rx-(BOX*round(delta_rx/BOX)));
#     delta_ry=(delta_ry-(BOX*round(delta_ry/BOX)));
#     edges_to_draw.append([(x1,y1),(x1-delta_rx,y1-delta_ry)])
#     colors.append(color)
#     
##color=(0.,1.,0.,1)
##for re in redundant_edges:
##    #print re
##    x1=CONFIG[re[0]][2]
##    y1=CONFIG[re[0]][3]
##    x2=CONFIG[re[1]][2]
##    y2=CONFIG[re[1]][3]
##    delta_rx=x1-x2
##    delta_ry=y1-y2
##
##    delta_rx=(delta_rx-(tilt*round(delta_ry/BOX)));
##    delta_rx=(delta_rx-(BOX*round(delta_rx/BOX)));
##    delta_ry=(delta_ry-(BOX*round(delta_ry/BOX)));
##    edges_to_draw.append([(x1,y1),(x1-delta_rx,y1-delta_ry)])
##    colors.append(color)
#lc = mc.LineCollection(edges_to_draw,colors=colors)
#ax.add_collection(lc)
#plt.savefig('network_'+(sys.argv[3])+'.pdf')
#plt.close(fig)
#########################################################################################################
###########################CODE TO DRAW CONFIGS##########################################################
################################END######################################################################
#########################################################################################################
#plt.show()
