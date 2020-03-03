import matplotlib.pyplot as plt

def view( coords, elt2verts, solGlob, colmap=plt.cm.Accent, title=None, visuMesh = True ) :
    fig = plt.figure()
    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.gca().set_aspect('equal')
    x = coords[:,0]
    y = coords[:,1]
    plt.tripcolor(x, y, elt2verts, solGlob, shading = 'gouraud', cmap=colmap)
    if visuMesh :
        plt.triplot(x, y, elt2verts, 'bo-', color='black')
    plt.axis('off')
    plt.colorbar()
    plt.show()
    
