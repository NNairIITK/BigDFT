#V_Sim resources file v3.0
#====================

#WARNING: this file format is DIFFERENT from that for
#standard v_sim version <= 2.x

#Line beginning with a # are not parsed.

#The only "useful" lines must have the following contents
#several two or more lines patterns:
#resource_name:
#values separeted by blank characters

#The following resource names are valid :
# atomic_radius_shape
# spin_resources
# pairWire_width
# pairWire_pairWidth
# pairWire_linkWidth
# pairWire_linkStipple
# cylinder_colorType
# pairCylinder_radius
# pairCylinder_pairRadius
# pairCylinder_linkRadius
# pairs_are_on
# pair_data
# pair_link
# pairs_favoriteMethod
# isosurfaces_drawIntra
# isosurface_property
# isosurface_color
# isosurface_properties
# nodeDisplacement_arrow
# nodeDisplacement_minThreshold
# nodeDisplacement_lblThreshold
# nodeDisplacement_factor
# highlight_radiusFactor
# scales_are_on
# scales_color
# scales_line_width
# scale_definition
# scales_line_stipple
# opengl_theta_phi_omega
# opengl_xs_ys
# opengl_gross
# opengl_d_red
# fog_is_on
# fog_color_is_specific
# fog_specific_color
# fog_start_end
# backgroundColor_color
# box_is_on
# box_color
# box_line_width
# box_line_stipple
# box_show_lengths
# axes_are_on
# axes_color
# axes_line_width
# axes_line_stipple
# legend_is_on
# material
# element_color
# element_is_rendered
# element_properties

# The radius of the element and its shape, a real > 0. & [Sphere Cube Elipsoid Point]
atomic_radius_shape:
    Si      1.175 Sphere
atomic_radius_shape:
    Al      1.430 Sphere
atomic_radius_shape:
    Ag      1.430 Sphere
atomic_radius_shape:
    Au      1.437 Sphere
atomic_radius_shape:
    Ge      1.225 Sphere
atomic_radius_shape:
    O      1.350 Sphere
atomic_radius_shape:
    C      0.600 Sphere
atomic_radius_shape:
    Cd      1.490 Sphere
atomic_radius_shape:
    Pd      1.375 Sphere
atomic_radius_shape:
    Te      1.430 Sphere
atomic_radius_shape:
    Co      1.250 Sphere
atomic_radius_shape:
    Fe      1.240 Sphere
atomic_radius_shape:
    Cu      1.250 Sphere
atomic_radius_shape:
    H      0.500 Sphere
atomic_radius_shape:
    Ni      1.243 Sphere
atomic_radius_shape:
    Pt      1.385 Sphere
atomic_radius_shape:
    g      0.250 Point
atomic_radius_shape:
    G      0.500 Point
atomic_radius_shape:
    N      0.392 Sphere

# Global or element resource for rendering spin module
spin_resources:
   spin_global_color_cone 0.000000 0.000000
spin_resources:
   spin_global_color_wheel 0.000000
spin_resources:
   spin_global_hiding_mode never
spin_resources:
   spin_global_atomic 0
spin_resources:
   spin_global_modulus 0

# This value is the width for all pairs drawn ; 0 < integer < 10
pairWire_width:
    2
# Widths detail for each drawn link ; 0 < integer < 10
pairWire_linkWidth:
    N N  2.100 2.200  6
pairWire_linkStipple:
    N N  2.100 2.200  61680

# It chooses the colors of the cylinders according differents criterion ; 0 <= integer < 2
cylinder_colorType:
    0
# This value is the default radius of the pairs drawn as cylinders ; 0 < real < 10
pairCylinder_radius:
    0.150000
# This value is the radius for specific pairs drawn as cylinders ; element1 elemen2 0 < real < 10

# Ask the opengl engine to draw pairs between elements ; boolean 0 or 1
pairs_are_on:
    1
# Favorite method used to render files ; chain ('Wire pairs', 'Cylinder pairs')
pairs_favoriteMethod:
    Wire pairs
# Draw a link between [ele1] [ele2] [0. <= dmin] [0. <= dmax]
#                     [0. <= RGB <= 1.]x3 [bool: drawn] [bool: printLength] [string: method]
pair_link:
    Au Ni 0.000 0.000
    1.000 0.600 0.200  1  0  Wire pairs
pair_link:
    Ni Ni 0.000 4.100
    0.750 0.400 0.200  1  0  Wire pairs
pair_link:
    Si Si 2.200 2.500
    0.750 0.400 0.200  1  0  Wire pairs
pair_link:
    Au Au 0.000 4.100
    1.000 0.600 0.200  1  0  Wire pairs
pair_link:
    N N 2.100 2.200
    1.000 1.000 1.000  1  0  Wire pairs

# Choose if the interior is drawn in color inverse ; a boolean (0 or 1)
isosurfaces_drawIntra:
    0

# Define the colour of one surface ; 4 floats (RGBA) 5 floats (material)
# Define some surface properties ; rendered (0 or 1) sensitive to planes (0 or 1)

# Describe the arrow to be drawn ; four floats (tail lg. tail rd. head lg. head rd., negative values means auto)
nodeDisplacement_arrow:
    -1.000000 -1.000000 -1.000000 -1.000000
# Choose the factor to draw arrows in geometry differences ; float (negative means auto)
nodeDisplacement_factor:
    -1.000000
# Choose the minimum value for drawn arrows in geometry differences ; float (ratio threshold if between -1 and 0)
nodeDisplacement_minThreshold:
    -0.100000
# Choose the minimum value for labels in geometry differences ; float (ratio threshold if between -1 and 0)
nodeDisplacement_lblThreshold:
    -0.900000
# Give the factor for the highlight radius ; one float (> 1.)
highlight_radiusFactor:
    1.250000
# Control if scales are drawn ; boolean (0 or 1)
scales_are_on:
    0
# Define the color RGBA of all scales ; four floating point values (0. <= v <= 1.)
scales_color:
    0.000 0.000 0.000
# Define the width of the lines of all scales ; one floating point value (1. <= v <= 10.)
scales_line_width:
       1
# Define the stipple pattern of the lines of all scales ; one integer value (0 <= v <= 65535)
scales_line_stipple:
    65535
# Define the position, the direction, the length and the legend of a scale ; position[3] direction[3] length legend
scale_definition:
    0 0 0  1 0 0  5  [auto]

# 2 real values (degrees) for user orientation with respect to sample
opengl_theta_phi_omega:
       68.500   -49.700     0.000
# 2 real values for image position with respect to [0.0, 1.0]x[0.0, 1.0] window
opengl_xs_ys:
        0.500     0.500
# gross factor (must be real > 0.0)
opengl_gross:
        4.177
# reduced perspective distance (must be real > 1.0)
opengl_d_red:
        5.000

# Control if the fog is used ; boolean (0 or 1)
fog_is_on:
    1
# Control if the fog uses a specific color ; boolean (0 or 1)
fog_color_is_specific:
    0
# Define the color of the fog ; four floating point values (0. <= v <= 1.)
fog_specific_color:
    0.000 0.000 0.000 0.000
# Define the position of the fog ; two floating point values (0. <= v <= 1.)
fog_start_end:
    0.447 0.631

# Set the background of the background ; four floating point values (0. <= v <= 1.)
backgroundColor_color:
    0.538 0.667 0.699 1.000

# Control if a box is drawn around the rendering area ; boolean (0 or 1)
box_is_on:
    0
# Define the color of the box ; three floating point values (0. <= v <= 1.)
box_color:
    0.823 0.984 0.849
# Define the width of the lines of the box ; one integer (1. <= v <= 10.)
box_line_width:
       2
# Dot scheme detail for the lines of the box (main and expanded) ; 0 < 2 integers < 2^16
box_line_stipple:
    65535 65280
# Print the box lengths ; boolean (0 or 1)
box_show_lengths:
    0

# Control if the axes are drawn ; boolean (0 or 1)
axes_are_on:
    0
# Define the color of the axes ; three floating point values (0. <= v <= 1.)
axes_color:
    0.823 0.884 0.849
# Define the width of the lines of the axes ; one floating point values (1. <= v <= 10.)
axes_line_width:
       1
# Dot scheme detail for the lines of the axes ; 0 < integer < 2^16
axes_line_stipple:
    65535

# Control if the legend is drawn ; boolean (0 or 1)
legend_is_on:
    1
# Codes the main color in RedGreenBlueAlpha formatand the light effects on material, nine floats between 0. and 1.
element_color:
    Si 0.000 1.000 0.200 1.000   0.20 0.54 0.69 0.50 0.20
element_color:
    Al 0.600 0.600 1.000 1.000   0.25 0.80 0.50 0.70 0.00
element_color:
    Ag 1.000 1.000 1.000 1.000   0.20 0.60 0.50 1.00 0.00
element_color:
    Au 1.000 0.750 0.040 1.000   0.54 0.54 0.16 0.52 0.00
element_color:
    Ge 0.000 1.000 0.800 1.000   0.20 0.54 0.69 0.50 0.20
element_color:
    O 1.000 0.200 0.200 1.000   0.60 0.60 0.00 0.00 0.00
element_color:
    C 0.300 0.300 0.300 1.000   0.20 0.54 0.69 0.50 0.20
element_color:
    Cd 0.000 0.800 1.000 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    Pd 1.000 1.000 0.800 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    Te 0.900 0.400 0.100 1.000   0.20 0.54 0.69 0.50 0.20
element_color:
    Co 1.000 0.500 0.800 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    Fe 1.000 0.500 0.000 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    Cu 1.000 0.560 0.340 1.000   0.20 0.70 0.10 0.24 0.15
element_color:
    H 1.000 1.000 1.000 1.000   0.60 0.60 0.00 0.00 0.00
element_color:
    Ni 1.000 0.300 0.300 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    Pt 1.000 0.800 0.800 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    g 1.000 0.750 0.040 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    G 0.600 0.600 1.000 1.000   0.50 0.00 0.12 0.50 0.20
element_color:
    N 0.000 1.000 0.200 1.000   0.25 0.25 0.25 0.25 0.25
# Define some properties ; rendered (0 or 1) masked(0 or 1).
element_properties:
    Si 1 1
element_properties:
    Ag 1 1
element_properties:
    Ge 1 1
element_properties:
    C 1 1
element_properties:
    Pd 1 1
element_properties:
    Co 1 1
element_properties:
    Cu 1 1
element_properties:
    Ni 1 1
element_properties:
    g 1 1
element_properties:
    N 1 1

