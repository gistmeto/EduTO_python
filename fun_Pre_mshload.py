import numpy as np

def fun_Pre_mshload(input):
    msh = {}
    flag_category = 0
    list_category = ["$MeshFormat","$EndMeshFormat","$PhysicalNames","$EndPhysicalNames","$Entities","$EndEntities","$Nodes","$EndNodes","$Elements","$EndElements"]
    flag_state = 0
    flag_debug_print = False

    with open(input["modelname"]+'.msh', 'r') as f:
        lines = f.readlines()
        lines = [line.rstrip('\n') for line in lines]   # remove \n
        
        for line in lines:        
            if line[0] == '$': ## Flag for parsing
                # print(line)
                for i in range(len(list_category)):
                    if line == list_category[i]:
                        if flag_debug_print:
                            print(line)
                            print(i)
                        flag_category = i
                        flag_state = 0
                        # print(flag_category)
                        
                        continue
            
            else: ## Data parsing
                if flag_category == 2: # Physical name
                    if flag_state == 0: # Initial state
                        if flag_debug_print:
                            print('Physical data parsing')
                        
                        if len(line) == 1:
                            num_phys_entities = int(line)
                            flag_state = 1

                            phys_tag_str = {}
                            phys_tag_str["RollerU"] = []
                            phys_tag_str["RollerV"] = []
                            phys_tag_str["RollerW"] = []
                            phys_tag_str["Fixed"] = []
                            phys_tag_str["Free"] = []
                            phys_tag_str["Traction"] = []
                            phys_tag_str["Body"] = []       # TODO
                            phys_tag_str["Design"] = []
                            phys_tag_str["NonDesign"] = []
                            phys_tag_str["External"] = []

                            bc = {}
                            bc["r"] = {}
                            bc["r"]["u"] = {}
                            bc["r"]["v"] = {}
                            bc["r"]["w"] = {}
                            bc["f"] = {}

                            bcru_ind = 0
                            bcrv_ind = 0
                            bcrw_ind = 0
                            bcf_ind = 0
                            
                            lc = {}
                            lc["t"] = {}
                            lc["f"] = {}
                            lc["b"] = {}

                            lcf_ind = 0
                            lct_ind = 0
                            lcb_ind = 0

                            ind_phys_entities = 0

                    elif flag_state == 1:
                        X = [idx for idx in line.split(' ')]

                        if X[2].startswith('"RollerU"'):
                            ind_phys_entities += 1
                            Boundary_roller_u_list = []
                            phys_tag_str["RollerU"].append(int(X[1]))                            
                        elif X[2].startswith('"RollerV"'):
                            ind_phys_entities += 1
                            Boundary_roller_v_list = []
                            phys_tag_str["RollerV"].append(int(X[1]))
                        elif X[2].startswith('"RollerW"'):
                            ind_phys_entities += 1
                            Boundary_roller_w_list = []
                            phys_tag_str["RollerW"].append(int(X[1]))
                        elif X[2].startswith('"Fixed"'):
                            ind_phys_entities += 1
                            Boundary_fixed_list = []
                            phys_tag_str["Fixed"].append(int(X[1]))
                        elif X[2].startswith('"Free"'):
                            ind_phys_entities += 1
                            Boundary_Free_list = []
                            phys_tag_str["Free"].append(int(X[1]))
                        elif X[2].startswith('"Traction"'):
                            ind_phys_entities += 1
                            Boundary_Traction_list = []
                            phys_tag_str["Traction"].append(int(X[1]))
                        elif X[2].startswith('"Body"'):
                            ind_phys_entities += 1
                            Body_list = []
                            phys_tag_str["Body"].append(int(X[1]))
                        elif X[2].startswith('"Design"'):
                            ind_phys_entities += 1
                            Design_list = []
                            phys_tag_str["Design"].append(int(X[1]))
                        elif X[2].startswith('"NonDesign"'):
                            ind_phys_entities += 1
                            NonDesign_list = []
                            phys_tag_str["NonDesign"].append(int(X[1]))
                        elif X[2].startswith('"External"'):
                            ind_phys_entities += 1
                            External_list = []
                            phys_tag_str["External"].append(int(X[1]))

                    if ind_phys_entities >= num_phys_entities:
                        flag_state = 0

                elif flag_category == 4: # Entities
                    if flag_state == 0:
                        X = [int(idx) for idx in line.split(' ')]

                        num_point = X[0]
                        num_line = X[1]
                        num_surf = X[2]
                        num_vol = X[3]
                        if 'phys_tag_str' in locals():
                            if (len(phys_tag_str) != 0):
                                flag_state = 1
                                flag_entity = 0
                                ind_entity = 0

                    elif flag_state == 1:
                        if flag_entity == 0: # Point
                            ind_entity += 1
                            if ind_entity >= num_point:
                                flag_entity = 1
                                ind_entity = 0

                        elif flag_entity == 1: # Line
                            ind_entity += 1
                            if ind_entity >= num_line:
                                flag_entity = 2
                                ind_entity = 0

                        elif flag_entity == 2: # Surface
                            ind_entity += 1

                            X = [idx for idx in line[:-1].split(' ')]

                            if int(X[7]) >= 1:
                                if int(X[8]) in phys_tag_str["RollerU"]:
                                    Boundary_roller_u_list.append(int(X[0]))
                                elif int(X[8]) in phys_tag_str["RollerV"]:
                                    Boundary_roller_v_list.append(int(X[0]))
                                elif int(X[8]) in phys_tag_str["RollerW"]:
                                    Boundary_roller_w_list.append(int(X[0]))
                                elif int(X[8]) in phys_tag_str["Fixed"]:
                                    Boundary_fixed_list.append(int(X[0]))
                                elif int(X[8]) in phys_tag_str["Free"]:
                                    Boundary_Free_list.append(int(X[0]))
                                elif int(X[8]) in phys_tag_str["Traction"]:
                                    Boundary_Traction_list.append(int(X[0]))

                            if ind_entity >= num_surf:
                                flag_entity = 3
                                ind_entity = 0

                        elif flag_entity == 3: # Volume
                            ind_entity += 1

                            X = [idx for idx in line[:-1].split(' ')]

                            for i in range(int(X[7])):
                                if int(X[8+i]) in phys_tag_str["Body"]:
                                    Body_list.append(int(X[0]))
                                elif int(X[8+i]) in phys_tag_str["Design"]:
                                    Design_list.append(int(X[0]))
                                elif int(X[8+i]) in phys_tag_str["NonDesign"]:
                                    NonDesign_list.append(int(X[0]))
                                elif int(X[8+i]) in phys_tag_str["External"]:
                                    External_list.append(int(X[0]))                            

                            if ind_entity >= num_surf:
                                flag_entity = 4
                                ind_entity = 0

                elif flag_category == 6: #Nodes
                    if flag_state == 0: # Initial state
                        if flag_debug_print:
                            print('Node data parsing')
                        # print(line)
                        X = [int(idx) for idx in line.split(' ')]
                        # print(X)
                        if len(X) != 4:
                            print("error")

                        num_node = X[1] # Total number of node
                        p = np.zeros((num_node,3))    # Node coordinates data
                        flag_state = 1

                    elif flag_state == 1:   # Info data parsing
                        if flag_debug_print:
                            print('Info data parsing')

                        X = [int(idx) for idx in line.split(' ')]
                        num_loop = X[3] # Number of node on this block
                        if flag_debug_print:
                            print(num_loop)
                            
                        if num_loop == 0:
                            continue
                        else:
                            flag_state = 2
                            ind_inner_loop = 0

                            ind_list = np.zeros((num_loop,1),dtype=int)

                    elif flag_state == 2:   # ind data parsing
                        if flag_debug_print:
                            print('ind data parsing')

                        ind_inner_loop += 1

                        ind_list[ind_inner_loop-1] = int(line)-1

                        if ind_inner_loop >= num_loop:  # If ind data parsing is end
                            flag_state = 3
                            ind_inner_loop = 0

                    elif flag_state == 3: # point data parsing
                        if flag_debug_print:
                            print('Point data parsing')

                        X = [float(idx) for idx in line.split(' ')]
                        ind_inner_loop = ind_inner_loop + 1
                        try:
                            p[ind_list[ind_inner_loop-1],:] = X[0:3]
                            if flag_debug_print:
                                print(p[ind_list[ind_inner_loop-1],:])
                        except:
                            print('except ind : %d' %(ind_inner_loop))
                        
                        if ind_inner_loop >= num_loop:
                            flag_state = 1
                            ind_inner_loop = 0
                            

                    elif flag_state == 99: # Error exception
                        ind_inner_loop += 1

                        if ind_inner_loop >= num_loop: # if ind data parsing is ended
                            flag_state = 1
                            ind_inner_loop = 0
                
                elif flag_category == 8: # $Elements
                    if flag_state == 0: # Initial state
                        if flag_debug_print:
                            print('Elements data parsing')

                        X = [int(idx) for idx in line.split(' ')]
                        if len(X) != 4:
                            print("error")
                        flag_state = 1

                        t = np.empty((0,5),dtype=int)
                        num_elem_list = []
                    
                    elif flag_state == 1:   # Info data parsing
                        # pri2nt(line)
                        try:
                            X = [int(idx) for idx in line.split(' ')]
                        except:
                            print(line)
                        num_loop = X[3]
                        ind_inner_loop = 0

                        flag_condition = 0

                        if X[0] == 2: # 2D element
                            if X[2] == 2:
                                if 'Boundary_roller_u_list' in locals():
                                    if X[1] in Boundary_roller_u_list:
                                        flag_state = 2
                                        flag_condition = 'U'
                                        bcru_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')  
                                if 'Boundary_roller_v_list' in locals():
                                    if X[1] in Boundary_roller_v_list:
                                        flag_state = 2
                                        flag_condition = 'V'
                                        bcrv_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')  
                                if 'Boundary_roller_w_list' in locals():
                                    if X[1] in Boundary_roller_w_list:
                                        flag_state = 2
                                        flag_condition = 'W'
                                        bcrw_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')    
                                if 'Boundary_fixed_list' in locals():
                                    if X[1] in Boundary_fixed_list:
                                        flag_state = 2
                                        flag_condition = 'F'
                                        bcf_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')  
                                if 'Boundary_Free_list' in locals():
                                    if X[1] in Boundary_Free_list:
                                        flag_state = 2
                                        flag_condition = 'N'
                                        lcf_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')   
                                if 'Boundary_Traction_list' in locals():
                                    if X[1] in Boundary_Traction_list:
                                        flag_state = 2
                                        flag_condition = 'T'
                                        lct_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,3),dtype='int')  
                                
                                    flag_state = 2
                                    num_elem = X[3]
                                    t_temp = np.zeros((num_elem,3),dtype='int')

                                if flag_condition == 0:
                                    flag_state = 98
                                    num_elem = X[3]


                        elif X[0] == 3: #
                            # pri2nt(line)
                            X = [int(idx) for idx in line.split(' ')]
                            num_loop = X[3]
                            ind_inner_loop = 0

                            flag_condition = 0
                            
                            if X[2] == 4: # 3D tetrahedron element
                                if 'Body_list' in locals():
                                    if X[1] in Body_list:
                                        flag_state = 3
                                        flag_condition = 'Body'
                                        lcb_ind += 1
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,4),dtype='int')  
                                if 'Design_list' in locals():
                                    if X[1] in Design_list:
                                        flag_state = 3
                                        flag_condition = 'Design'
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,5),dtype='int')  
                                if 'NonDesign_list' in locals():
                                    if X[1] in NonDesign_list:
                                        flag_state = 3
                                        flag_condition = 'NonDesign'
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,5),dtype='int')    
                                if 'External_list' in locals():
                                    if X[1] in External_list:
                                        flag_state = 3
                                        flag_condition = 'External'
                                        num_elem = X[3]
                                        t_temp = np.zeros((num_elem,5),dtype='int')                          

                        else:
                            flag_state = 98
                            num_elem = X[3]

                    elif flag_state == 2:   # 
                        ind_inner_loop += 1
                        X = [int(idx) for idx in line[:-1].split(' ')] # line[:-1] to remove last blank
                        try:
                            t_temp[ind_inner_loop-1,:] = X[1:]
                        except:
                            print(X)


                        if ind_inner_loop >= num_loop:
                            flag_state = 1
                            ind_inner_loop = 0

                            if flag_condition == 'U':
                                bc['r']['u'][bcru_ind] = t_temp
                            elif flag_condition == 'V':
                                bc['r']['v'][bcrv_ind] = t_temp
                            elif flag_condition == 'W':
                                bc['r']['w'][bcrw_ind] = t_temp
                            elif flag_condition == 'F':
                                bc['f'][bcf_ind] = t_temp
                            elif flag_condition == 'N':
                                lc['f'][lcf_ind] = t_temp
                            elif flag_condition == 'T':
                                lc['t'][lct_ind] = t_temp
                                
                    elif flag_state == 3: # 3D Constan Tetrahedron element conectivity data parsing
                        ind_inner_loop += 1
                        X = [int(idx) for idx in line[:-1].split(' ')]

                        try: # TODO : Current version is only consier tetrahydral element
                            t_temp[ind_inner_loop-1,0:4] = X[1:]     # 1 2 3 4
                            
                            if flag_condition == "Design":
                                t_temp[ind_inner_loop-1,-1] = 1          # Material property
                            elif flag_condition == "NonDesign":
                                t_temp[ind_inner_loop-1,-1] = 2          # Material property
                            elif flag_condition == "External":
                                t_temp[ind_inner_loop-1,-1] = 3          # Material property
                        except:
                            print('except ind : %d' %(ind_inner_loop))

                        if ind_inner_loop >= num_elem:
                            flag_state = 1
                            ind_inner_loop = 0

                            if flag_condition != 'Body':
                                t = np.vstack((t,t_temp))
                            
                    elif flag_state == 98:   # 
                        ind_inner_loop = ind_inner_loop + 1
                        X = [int(idx) for idx in line[:-1].split(' ')]

                        if ind_inner_loop >= num_elem:
                            flag_state = 1
                            ind_inner_loop = 0

                    elif flag_state == 99:
                        ind_inner_loop = ind_inner_loop + 1

                        if ind_inner_loop >= num_loop:
                            flag_state = 1
                            ind_inner_loop = 0

    msh = {"X" : p, "IX" : t, "bc" : bc, "lc" : lc}
    msh["ne"] = msh["IX"].shape[0]
    
    A = np.zeros((msh["ne"],))
    nx = np.transpose(np.vstack([msh["X"][msh["IX"][:,0]-1,[0]],msh["X"][msh["IX"][:,1]-1,[0]],msh["X"][msh["IX"][:,2]-1,[0]]]))
    ny = np.transpose(np.vstack([msh["X"][msh["IX"][:,0]-1,[1]],msh["X"][msh["IX"][:,1]-1,[1]],msh["X"][msh["IX"][:,2]-1,[1]]]))

    for e in range(msh["ne"]):
        x = nx[e]
        y = ny[e]
        Ae = 1./2 * np.linalg.det(np.vstack(([np.ones((3,)),x,y]))) #TODO

        A[e] = np.abs(Ae)

    msh["WA"] = A

    L = np.zeros(lc["t"][1].shape[0])
    nx = np.transpose(np.vstack([msh["X"][lc["t"][1][:,0]-1,[0]],msh["X"][lc["t"][1][:,1]-1,[0]]]))
    ny = np.transpose(np.vstack([msh["X"][lc["t"][1][:,0]-1,[1]],msh["X"][lc["t"][1][:,1]-1,[1]]]))

    for e in range(lc["t"][1].shape[0]):
        x = nx[e,:]
        y = ny[e,:]
        Le = np.sqrt((x[0]-x[1])**2+(y[0]-y[1])**2)

        L[e] = Le

    msh["WL"] = L

    return msh