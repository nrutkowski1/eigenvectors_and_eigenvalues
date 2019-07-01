def cubic(a, b, c, d):
    
    p = (b**2 - 3.*a*c)
    q = (9.*a*b*c - 2*(b**3) - 27*(a**2)*d)
    n = (27*(p**3) / (q**2))
    theta = np.arccos(q / ((2.*p)*np.sqrt(p)))
    root1 = (-b + 2.*np.cos(theta / 3.)*np.sqrt(p)) / (3.*a)
    root2 = (-b + 2.*np.cos((theta / 3.) + (120. * (np.pi / 180)))*np.sqrt(p)) / (3.*a)
    root3 = (-b + 2.*np.cos((theta / 3.) + (240. * (np.pi / 180)))*np.sqrt(p)) / (3.*a)
    return [root1, root2, root3]       

def eigvals(matrix_3x3):
    
    vals = np.zeros(3)
    a,b,c,d,e,f,g,h,i = matrix_3x3.flatten()
    tr_m3x3 = a + e + i
    matrix_3x3_2 = np.dot(matrix_3x3, matrix_3x3)
    tr2_m3x3 =  matrix_3x3_2[0][0] + matrix_3x3_2[1][1] + matrix_3x3_2[2][2]
    det_m3x3 = (a*e*i - a*f*h) + (b*f*g - b*d*i) + (c*d*h - c*e*g)
    vals = cubic(-1, tr_m3x3, -0.5*(tr_m3x3**2 - tr2_m3x3), det_m3x3)
    vals = sorted(vals, key=abs, reverse=True)
    return vals
    
    def quadratic(a,b,c):
    sqrt_part = np.lib.scimath.sqrt(b**2 - 4*a*c)
    root1 = (-b + sqrt_part) / (2 * a)
    root2 = (-b - sqrt_part) / (2 * a)
    return root1, root2

def eigvals2x2(matrix_2x2):
    vals = np.zeros(2)
    a,b,c,d = matrix_2x2.flatten()
    vals[:] = quadratic(1.0, -(a+d), (a*d-b*c))
    vals = sorted(vals, key=abs, reverse=True)
    return vals

def eigvecs(matrix, vals):
    
    vecs = []
    
    for v in vals:
        
        
        a, b, c, d, e, f, g, h, i = np.array(matrix).flatten()
        row0 = [a - v, b, c]
        row1 = [d, e - v, f]
        row2 = [g, h, i - v]
        
        r0xr1 = [row0[1] * row1[2] - row0[2] * row1[1], 
                 row0[2] * row1[0] - row0[0] * row1[2],
                 row0[0] * row1[1] - row0[1] * row1[0]]
        r0xr2 = [row0[1] * row2[2] - row0[2] * row2[1], 
                 row0[2] * row2[0] - row0[0] * row2[2],
                 row0[0] * row2[1] - row0[1] * row2[0]]
        r1xr2 = [row1[1] * row2[2] - row1[2] * row2[1], 
                 row1[2] * row2[0] - row1[0] * row2[2],
                 row1[0] * row2[1] - row1[1] * row2[0]]
        d0 = r0xr1[0]*r0xr1[0] + r0xr1[1]*r0xr1[1] + r0xr1[2]*r0xr1[2]
        d1 = r0xr2[0]*r0xr2[0] + r0xr2[1]*r0xr2[1] + r0xr2[2]*r0xr2[2]
        d2 = r1xr2[0]*r1xr2[0] + r1xr2[1]*r1xr2[1] + r1xr2[2]*r1xr2[2]
        imax = 0
        dmax = d0
        
        if d1 > dmax:
            dmax = d1
            imax = 1
        if d2 > dmax:
            imax = 2
        if imax == 0:
            vecs.append(r0xr1 / np.sqrt(d0))
        if imax == 1:
            vecs.append(r0xr2 / np.sqrt(d1))
        if imax == 2:
            vecs.append(r1xr2 / np.sqrt(d2))
            
    return vecs
