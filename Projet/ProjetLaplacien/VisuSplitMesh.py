import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

def view( coords, elt2verts, nbDoms, partition,
          indInterfNodes = None, visuIndElts = True,
          visuIndVerts = True, colmap = plt.cm.Accent, title = None,
          isPartVert = True) :
    fntcol = '#004000'
    quot = max(1, nbDoms-1)
    fig = plt.figure()
    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.gca().set_aspect('equal')
    x = coords[:,0]
    y = coords[:,1]
    if (partition.shape[0] == coords.shape[0]) and isPartVert :
        # C'est une partition par noeud :
        plt.tripcolor(x, y, elt2verts, partition, shading = 'gouraud', cmap=colmap)
        plt.triplot(x, y, elt2verts, 'bo-', color='black')
    else :
        # C'est une partition par element :
        isPartVert = False
        plt.tripcolor(x, y, elt2verts, facecolors = partition,
                      edgecolors = 'k', cmap=colmap)
    if visuIndElts :
        for ie in range(elt2verts.shape[0]) :
            e = elt2verts[ie,:]
            bary = [(x[e[0]]+x[e[1]]+x[e[2]])/3., (y[e[0]]+y[e[1]]+y[e[2]])/3.]
            plt.text(bary[0],bary[1], u"%i"%ie, fontsize=10, fontweight='light',
                     style='italic', ha='center', va='center', color=fntcol)
    if visuIndVerts :
        fntcol = 'black'
        iv = 0
        if isPartVert :
            for v in coords :
                szfnt = 11
                wfont = 'light'
                plt.text(v[0],v[1], u"%i"%iv,ha='center', va='center', 
                         fontsize=szfnt, color=fntcol, fontweight=wfont,
                         backgroundcolor=colmap(partition[iv]/quot))
                iv += 1
        else :
            szfnt = 10
            wfont = 'light'
            mask = np.ones((coords.shape[0],), np.short)
            ie = 0
            for e in elt2verts :
                d = partition[ie]
                bcolor = colmap(d/quot)            
                if mask[e[0]] == 1 :
                    plt.text(x[e[0]],y[e[0]], u"%i"%e[0],ha='center', va='center',
                             fontsize=szfnt, color=fntcol, fontweight=wfont, backgroundcolor=bcolor)
                    mask[e[0]] = 0
                if mask[e[1]] == 1 :
                    plt.text(x[e[1]],y[e[1]], u"%i"%e[1],ha='center', va='center',
                             fontsize=szfnt, color=fntcol, fontweight=wfont, backgroundcolor=bcolor)
                    mask[e[1]] = 0
                if mask[e[2]] == 1 :
                    plt.text(x[e[2]],y[e[2]], u"%i"%e[2],ha='center', va='center',
                             fontsize=szfnt, color=fntcol, fontweight=wfont, backgroundcolor=bcolor)
                    mask[e[2]] = 0
                ie += 1

    if (indInterfNodes is not None) and visuIndVerts :
        for iv in indInterfNodes :
            v = coords[iv,:]
            szfnt = 12
            wfont = 'bold'
            d = partition[iv]
            bcolor = colmap(d/quot)            
            if not isPartVert :
                bcolor = 'white'
                fntcol = 'black'
            plt.text(v[0],v[1], u"%i"%iv,ha='center', va='center',
                     fontsize=szfnt, color=fntcol, fontweight=wfont,
                     backgroundcolor=bcolor)
    plt.axis('off')
    plt.show()
