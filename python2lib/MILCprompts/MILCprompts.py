# Procedures for generating MILC prompts
# Mostly by J. Simone.  Some additions by C. DeTar

import sys
import textwrap
from Cheetah.Template import Template

def base36(n):
    """Convert a positive integer to a base36 string."""
    alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    base = 36
    if n <= 0:
        return str(0)
    nb = ''
    while n != 0:
        n, i = divmod(n,base)
        nb = alphabet[i] + nb
        pass
    return nb

def oppMom(mom):
    """Take negative of momentum coordinates."""
    return [-p for p in mom]

class Geometry:
    """Lattice geometry, machine layout (SciDAC) and seed for random nuber generator."""
    _Template = """
    #== ${_classType} ==
    prompt ${prompt}
    nx ${dim[0]}
    ny ${dim[1]}
    nz ${dim[2]}
    nt ${dim[3]}
    #if $layout is not None:
    node_geometry #echo ' '.join(map(str,$layout.node))#
    ionode_geometry #echo ' '.join(map(str,$layout.io))# 
    #end if
    iseed ${seed}
    job_id ${jobID}"""
    def __init__(self,dim,seed,jobID,prompt,layout):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.dim = dim
        self.seed = seed
        self.jobID = jobID
        self.layout = layout
        self.prompt = prompt
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class Gauge:
    """The SU3 gauge field, gauge fixing and APE link smearing."""
    _Template = """
    #== ${_classType} ==
    #echo ' '.join($load)#
    u0 ${u0}
    ${gFix}
    #echo ' '.join($save)#
    staple_weight ${fatLink.weight}
    ape_iter ${fatLink.iter}
    coordinate_origin #echo ' '.join(map(str,$origin))#
    time_bc ${bc}"""
    def __init__(self,load,u0,gFix,save,fatLink,origin,bc):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.load = load
        self.u0 = u0
        self.gFix = gFix
        self.save = save
        self.fatLink = fatLink
        self.origin = origin
        self.bc = bc
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class PointSource:
    """Point base source."""
    _Template = """
    #== source ${id}: ${_classType} ==
    point
    field_type ${field_type}
    subset ${subset}
    origin #echo ' '.join(map(str,$origin))#
    #if $scaleFactor is not None:
    scale_factor ${scaleFactor}
    #end if
    source_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,origin,field_type,subset,scaleFactor,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.origin = origin
        self.field_type = field_type
        self.subset = subset
        self.scaleFactor = scaleFactor
        self.label = label
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        pass
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class RandomColorWallSource:
    """Random color wall base source."""
    _Template = """
    #== source ${id}: ${_classType} ==
    random_color_wall
    field_type ${field_type}
    subset ${subset}
    t0 ${tsrc}
    ncolor ${ncolor}
    momentum #echo ' '.join(map(str,$momentum))#
    #if $scaleFactor is not None:
    scale_factor ${scaleFactor}
    #end if
    source_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,tsrc,ncolor,field_type,subset,momentum,scaleFactor,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.tsrc = tsrc
        self.ncolor = ncolor
        self.field_type = field_type
        self.subset = subset
        self.momentum = momentum
        self.scaleFactor = scaleFactor
        self.label = label
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class VectorFieldSource:
    """Color vector field base source input from a file."""
    _Template = """
    #== source ${id}: ${_classType} ==
    vector_field
    field_type ${field_type}
    subset ${subset}
    origin #echo ' '.join(map(str,$origin))#
    #echo ' '.join($load)#
    ncolor ${ncolor}
    momentum #echo ' '.join(map(str,$momentum))#
    #if $scaleFactor is not None:
    scale_factor ${scaleFactor}
    #end if
    source_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,load,origin,ncolor,field_type,subset,momentum,scaleFactor,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.load = load
        self.origin = origin
        self.ncolor = ncolor
        self.field_type = field_type
        self.subset = subset
        self.momentum = momentum
        self.scaleFactor = scaleFactor
        self.label = label
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class DiracFieldSource:
    """Dirac (spin and color) base source."""
    _Template = """
    #== source ${id}: ${_classType} ==
    dirac_field
    field_type clover
    subset ${subset}
    origin #echo ' '.join(map(str,$origin))#
    #echo ' '.join($load)#
    nsource #echo 4 * ${ncolor}#
    momentum #echo ' '.join(map(str,$momentum))#
    #if $scaleFactor is not None:
    scale_factor ${scaleFactor}
    #end if
    source_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,load,origin,ncolor,subset,momentum,scaleFactor,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.load = load
        self.origin = origin
        self.ncolor = ncolor
        self.subset = subset
        self.momentum = momentum
        self.scaleFactor = scaleFactor
        self.label = label
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class RadialWavefunction:
    """Radial wavefunction source from a file. Either as a base or a derived source type.""" 
    _Template = """
    #== source ${id}: ${_classType} ==
    #if $startSource is not None:
    source ${startSource.id}
    #end if
    wavefunction
    #echo ' '.join($load)#
    #if $stride is not None:
    stride ${stride}
    #end if
    a ${afm}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,label,afm,stride,load,save,startSource=None):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.label = label
        self.afm = afm
        self.stride = stride
        self.load = load
        self.save = save
        self.startSource = startSource
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        if self.startSource is not None:
            depends.append(self.startSource._objectID)
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class FermilabRotation:
    """Apply a Fermilab rotation to a source."""
    _Template = """
    #== source ${id}: ${_classType} ==
    #if $startSource is not None:
    source ${startSource.id}
    #end if
    rotate_3D
    d1 ${d1}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,label,d1,save,startSource=None):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.label = label
        self.d1 = d1
        self.save = save
        self.startSource = startSource
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        if self.startSource is not None:
            depends.append(self.startSource._objectID)
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class FatCovariantGaussian:
    """Covariant Gaussian source from APE smeared links. Either as a base or a derived source type."""
    _Template = """
    #== source ${id}: ${_classType} ==
    #if $startSource is not None:
    source ${startSource.id}
    #end if
    fat_covariant_gaussian
    stride ${gparams.stride}
    r0 ${gparams.r0}
    source_iters ${gparams.iters}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,gparams,label,save,startSource=None):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.gparams = gparams
        self.label = label
        self.save = save
        self.startSource = startSource
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        if self.startSource is not None:
            depends.append(self.startSource._objectID)
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class Eigen:
    """Read eigenpairs from file."""
    _Template = """
    #== ${_classType} ==
    max_number_of_eigenpairs ${Nvecs}
    #if $Nvecs > 0:
    #echo ' '.join($load)#
    #echo ' '.join($save)#
    #end if"""
    def __init__(self,load, Nvecs, save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.Nvecs = Nvecs
        self.load = load
        self.save =save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        depends.append(self.source._objectID)
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass
    

class SolveKS:
    """su3_clov KS propagator solve."""
    _Template = """
    #== propagator ${id}: ${_classType} ==
    propagator_type KS
    mass ${mass}
    #if $naik_epsilon is not None:
    naik_term_epsilon ${naik_epsilon}
    #end if
    check ${check}
    error_for_propagator ${residual.L2}
    rel_error_for_propagator ${residual.R2}
    precision ${precision}
    momentum_twist #echo ' '.join(map(str,$twist))#
    source ${source.id}
    #echo ' '.join($load)#
    #echo ' '.join($save)#"""
    def __init__(self,mass,naik_epsilon,source,twist,load,save,residual,precision,check):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.mass = mass
        self.naik_epsilon = naik_epsilon
        self.source = source
        self.twist = twist
        self.load = load
        self.save = save
        self.residual = residual
        self.precision = precision
        self.check = check
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        depends.append(self.source._objectID)
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class SolveClover:
    """su3_clov Clover propagator solve."""
    _Template = """
    #== propagator ${id}: ${_classType} ==
    propagator_type clover
    kappa ${kappa}
    clov_c ${cSW}
    check ${check}
    error_for_propagator ${residual.L2}
    rel_error_for_propagator ${residual.R2}
    precision ${precision}
    momentum_twist #echo ' '.join(map(str,$twist))#
    source ${source.id}
    #echo ' '.join($load)#
    #echo ' '.join($save)#"""
    def __init__(self,kappa,cSW,source,twist,load,save,residual,precision,check):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.kappa = kappa
        self.cSW = cSW
        self.source = source
        self.twist = twist
        self.load = load
        self.save =save
        self.residual = residual
        self.precision = precision
        self.check = check
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        depends.append(self.source._objectID)
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class QuarkIdentitySink:
    """NOP fermion sink smearing."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    identity
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class RadialWavefunctionSink:
    """Smear fermion sink with a radial wavefunction input from file."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    wavefunction
    #echo ' '.join($load)#
    #if $stride is not None:
    stride ${stride}
    #end if
    a ${afm}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,label,afm,stride,load,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.afm = afm
        self.stride = stride
        self.load = load
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class FermilabRotateSink:
    """Apply a Fermilab rotation to a propagator."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    rotate_3D
    d1 ${d1}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,label,d1,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.d1 = d1
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class DiracExtSrcSink:
    """Restrict Dirac fermion propagator to time slice."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    ext_src_dirac
    gamma ${gamma}
    momentum #echo ' '.join(map(str,$momentum))#
    t0 ${time_slice}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,gamma,momentum,time_slice,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.gamma = gamma
        self.momentum = momentum
        self.time_slice = time_slice
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class KSExtSrcSink:
    """Restrict KS fermion propagator to time slice."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    ext_src_ks
    spin_taste_extend ${spin_taste_op}
    momentum #echo ' '.join(map(str,$momentum))#
    t0 ${time_slice}
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,spin_taste_op,momentum,time_slice,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.spin_taste_op = spin_taste_op
        self.momentum = momentum
        self.time_slice = time_slice
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class KSInverseSink:
    """Extended KS inverse sink operator."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    ks_inverse
    mass ${mass}
    #if $naik_epsilon is not None:
    naik_term_epsilon ${naik_epsilon}
    #end if
    u0 ${u0}
    max_cg_iterations ${maxCG.iters}
    max_cg_restarts ${maxCG.restarts}
    deflate ${deflate}
    error_for_propagator ${residual.L2}
    rel_error_for_propagator ${residual.R2}
    precision ${precision}
    momentum_twist #echo ' '.join(map(str,$twist))#
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,mass,naik_epsilon,u0,maxCG,deflate,residual,precision,twist,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.mass = mass
        self.naik_epsilon = naik_epsilon
        self.u0 = u0
        self.maxCG = maxCG
        self.deflate = deflate
        self.residual = residual
        self.precision = precision
        self.twist = twist
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class DiracInverseSink:
    """Extended Dirac inverse sink operator."""
    _Template = """
    #== quark ${id}: ${_classType} ==
    ${prop.type} ${prop.id}
    dirac_inverse
    kappa ${kappa}
    clov_c ${clov_c}
    u0 ${u0}
    max_cg_iterations ${maxCG.iters}
    max_cg_restarts ${maxCG.restarts}
    error_for_propagator ${residual.L2}
    rel_error_for_propagator ${residual.R2}
    precision ${precision}
    momentum_twist #echo ' '.join(map(str,$twist))#
    op_label ${label}
    #echo ' '.join($save)#"""
    def __init__(self,prop,kappa,clov_c,u0,maxCG,residual,precision,twist,label,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prop = prop
        self.label = label
        self.kappa = kappa
        self.clov_c = clov_c
        self.u0 = u0
        self.maxCG = maxCG
        self.residual = residual
        self.precision = precision
        self.twist = twist
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.prop._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class MesonSpectrum:
    """Meson spectrum specification."""
    _Template = """
    #== ${_classType} ==
    pair ${antiQuark.id} ${quark.id}
    spectrum_request meson
    #echo ' '.join($save)#
    r_offset #echo ' '.join(map(str,$relOffset))#
    number_of_correlators #echo len($npts)"""
    def __init__(self,antiQuark,quark,relOffset,npts,save):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.antiQuark = antiQuark
        self.quark = quark
        self.relOffset = relOffset
        self.npts = npts
        self.save = save
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        for c in self.npts:
            c.generate(ostream)
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.antiQuark._objectID)
        depends.append(self.quark._objectID)
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class MesonNpt:
    """A meson n-point function."""
    _Template = """correlator ${prefix} ${postfix} #echo ' '.join(map(str,$norm))# #echo ' '.join($gamma)# #echo ' '.join(map(str,$momentum))# #echo ' '.join($parity)"""
    def __init__(self,prefix,postfix,norm,gamma,momentum,parity):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.prefix = prefix
        self.postfix = postfix
        self.norm = norm
        self.gamma = gamma
        self.momentum = momentum
        self.parity = parity
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return

class _Cycle_su3_clov:
    def __init__(self):
        self.gauge = None
        self.maxIters = None
        self.maxRestarts = None
        self.bsource = list()
        self.msource = list()
        self.propagator = list()
        self.quark = list()
        self.spectrum = list()
        return
    def generate(self,ostream):
        self.gauge.generate(ostream)
        print>>ostream, 'max_cg_iterations', self.maxIters
        print>>ostream, 'max_cg_restarts', self.maxRestarts
        print>>ostream
        print>>ostream, '########### base sources ###############'
        print>>ostream
        print>>ostream, 'number_of_base_sources', len(self.bsource)
        for x in self.bsource:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '########### modified sources ###############'
        print>>ostream
        print>>ostream, 'number_of_modified_sources', len(self.msource)
        for x in self.msource:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '########### propagators ###############'
        print>>ostream
        print>>ostream, 'number_of_propagators', len(self.propagator)
        for x in self.propagator:
            x.generate(ostream)
        print>>ostream
        print>>ostream, '########### quarks ###############'
        print>>ostream
        print>>ostream, 'number_of_quarks', len(self.quark)
        for x in self.quark:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '########### mesons ###############'
        print>>ostream
        print>>ostream, 'number_of_pairings', len(self.spectrum)
        for x in self.spectrum:
            x.generate(ostream)
        return
    def dataflow(self):
        dinfo = list()
        dinfo.append(self.gauge.dataflow())
        for x in self.bsource:
            dinfo.append(x.dataflow())
            pass
        for x in self.msource:
            dinfo.append(x.dataflow())
            pass
        for x in self.propagator:
            dinfo.append(x.dataflow())
            pass
        for x in self.quark:
            dinfo.append(x.dataflow())
            pass
        for x in self.spectrum:
            dinfo.append(x.dataflow())
            pass
        return dinfo
    pass

class su3_clov:
    """
    The MILC/su3_clov application.
    
    :Parameters:
      - `participantName`: A unique label for the workflow.
      - `dim`: List of lattice dimensions, [nX, nY, nZ, nT].
      - `seed`: Random number generator seed, use seed='None' to skip.
      - `jobID`: A unique identifier string used for data provenance, usually set to PBS_JOBID.
      - `layout`: The SciDAC node and io layout = { node: [Lx, Ly, Lz, Lt], io: [Ix, Iy, Iz, It] }. Use layout = 'None' to skip.
      - `prompt`: Echo MILC prompts to stdout (=0), no echo (=1) or check input validity (=2).
    
    The 'su3_clov' application is a pipline:

    - initialize lattice size, rng.seed
    - loop:
        - initialize gauge field
        - gauge fix
        - fat link APE smear
        - define base sources
        - define derived sources
        - define propagator solves
        - define quarks (propagator sink treatments)
        - spectroscopy: n-point tie-ups

    """
    def __init__(self,participantName,dim,seed,jobID,layout=None,prompt=0):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.pName = participantName
        self.geometry = Geometry(dim,seed,jobID,prompt,layout)
        self.requires = list()
        self.produces = list()
        self.cycle = [ ]
        self.cycnt = -1
        return
    def newGauge(self,uspec):
        """Add the gauge definition."""
        self.cycle.append(_Cycle_su3_clov())
        self.cycnt += 1
        self.cycle[self.cycnt].gauge = uspec
        return uspec
    def setCGparams(self,maxIters,maxRestarts):
        """Specify CG maximum restarts and iterations between restarts."""
        self.cycle[self.cycnt].maxIters = maxIters
        self.cycle[self.cycnt].maxRestarts = maxRestarts
        return
    def addBaseSource(self,src):
        """Add a base source specification"""
        self.cycle[self.cycnt].bsource.append(src)
        return src
    def addModSource(self,src):
        """Add a source derived from one of the base sources"""
        self.cycle[self.cycnt].msource.append(src)
        return src
    def addProp(self,prop):
        """Add a propagator solve."""
        prop.type = 'propagator'
        self.cycle[self.cycnt].propagator.append(prop)
        return prop
    def addQuark(self,quark):
        """Add a sink treatment to a propagator."""
        quark.type = 'quark'
        self.cycle[self.cycnt].quark.append(quark)
        return quark
    def addSpectrum(self,spect):
        """Add a spectrum specification that depends upon pairs of quarks."""
        self.cycle[self.cycnt].spectrum.append(spect)
        return spect
    def _bind_indices(self):
        """Bind indices to sources, propagators and quarks."""
        for cy in self.cycle:
            nbs = len(cy.bsource)
            for idx, s in zip(range(nbs),cy.bsource):
                s.id = idx
                pass
            nms = len(cy.msource)
            for idx, s in zip(range(nbs,nbs+nms),cy.msource):
                s.id = idx
                pass
            for idx, p in zip(range(len(cy.propagator)),cy.propagator):
                p.id = idx
                pass
            for idx, q in zip(range(len(cy.quark)),cy.quark):
                q.id = idx
                pass
        return
    def generate(self,ostream=sys.stdout):
        """Write MILC prompts to output stream 'ostream'."""
        self._bind_indices()
        self.geometry.generate(ostream)
        for x in self.cycle:
            x.generate(ostream)
            pass
        return
    def dataflow(self):
        self._bind_indices()
        w = list()
        for x in self.cycle:
            w.append(x.dataflow())
            pass
        dinfo = { 'title': self.pName, 'workflow': w }
        return dinfo
    pass

#----

class KSsolveSet:
    """A set of KS solves that have a common source specification, momentum twist, and precision."""
    _Template = """
    #== ${_classType} ==
    set_type ${set_type}
    inv_type ${inv_type}
    max_cg_iterations ${maxCG.iters}
    max_cg_restarts ${maxCG.restarts}
    check ${check}
    momentum_twist #echo ' '.join(map(str,$twist))#
    precision ${precision}
    source ${source.id}"""
    def __init__(self,source,twist,check,set_type,inv_type,maxCG,precision):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.source = source
        self.twist = twist
        self.check = check
        self.set_type = set_type
        self.inv_type = inv_type
        self.maxCG = maxCG
        self.precision = precision
        self.propagator = list()
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    # container behavior
    def __len__(self):
        return len(self.propagator)
    def __iter__(self):
        return self.propagator.__iter__()
    def __add__(self,other):
        return [ x for x in self.propagator+other.propagator ]
    def addPropagator(self,prop):
        """Add a KSsolveElement object to the set of solves."""
        prop.parent = self
        prop.type = 'propagator'
        self.propagator.append(prop)
        return prop
    def generate(self,ostream):
        print>>ostream, self._template
        print>>ostream, 'number_of_propagators', len(self.propagator)
        for p in self.propagator:
            p.generate(ostream)
            pass
        return
    def dataflow(self):
        df = list()
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.source._objectID)
        df.append( { 'name': wname, 'depends': depends,
                     'requires': requires, 'produces': produces } )
        for p in self.propagator:
            df.append(p.dataflow())
            pass
        return df
    pass

class KSsolveElement:
    """Specification of a single KS solve in a KSsolveSet."""
    _Template = """
    #== propagator ${id}: ${_classType} ==
    mass ${mass}
    #if $naik is not None
    naik_term_epsilon ${naik}
    #end if
    #if $deflate is not None:
    deflate ${deflate}
    #end if
    error_for_propagator ${residual.L2}
    rel_error_for_propagator ${residual.R2}
    #echo ' '.join($load)#
    #echo ' '.join($save)#
    """
    def __init__(self,mass,naik,load,save,deflate,residual):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.id = None
        self.parent = None # the KSsolveSet
        self.mass = mass
        self.naik = naik
        self.load = load
        self.save = save
        self.deflate = deflate
        self.residual = residual
        self._template = Template(source=textwrap.dedent(self._Template),searchList=vars(self))
        return
    def generate(self,ostream):
        print>>ostream, self._template
        return
    def dataflow(self):
        wname = self._objectID
        depends = list()
        requires = list()
        produces = list()
        depends.append(self.parent._objectID)
        if len(self.load) > 1:
            requires.append(self.load[1])
            pass
        if len(self.save) > 1:
            produces.append(self.save[1])
            pass
        return { 'name': wname, 'depends': depends,
                 'requires': requires, 'produces': produces }
    pass

class ks_spectrum:
    """
    MILC application MILC/ks_spectrum application.
    The milc applications are designed as pipelines.
    The ks_spectrum pipline:
    
    - initialize lattice size, rng
    - loop:
        - init gauge field
        -  gauge fix
        -  link smear
        -  eigenpairs??
        -  pbp masses??
        -  define sources
        -  define KS propsets
        -  define quarks (sink treatments)
        -  meson spectroscopy
        -  baryon spectroscopy
        
    """
    def __init__(self,participantName,dim,seed,jobID,layout=None,prompt=0):
        self._classType = self.__class__.__name__
        self._objectID = self._classType+'_'+base36(id(self))
        self.pName = participantName
        self.geometry = Geometry(dim,seed,jobID,prompt,layout)
        self.requires = list()
        self.produces = list()
        self.cycle = [ ]
        self.cycnt = -1
        return
    def newGauge(self,uspec):
        """Add a gauge field specification."""
        self.cycle.append(_Cycle_ks_spectrum())
        self.cycnt += 1
        self.cycle[self.cycnt].gauge = uspec
        return uspec
    def newEigen(self,eigspec):
        """Add an eigensolve specification."""
        self.cycle[self.cycnt].eigen = eigspec
        return eigspec
    def addBaseSource(self,src):
        """Add a base source specification."""
        self.cycle[self.cycnt].bsource.append(src)
        return src
    def addModSource(self,src):
        """Add a modified source dependent upon one of the base sources."""
        self.cycle[self.cycnt].msource.append(src)
        return src
    def addPropSet(self,pset):
        """Add a KSsolveSet object."""
        self.cycle[self.cycnt].pset.append(pset)
        return pset
    def addQuark(self,quark):
        """Add a propagator sink treatment (a quark) to the workflow. The quark depends upon an KSsolveElement object."""
        quark.type = 'quark'
        self.cycle[self.cycnt].quark.append(quark)
        return quark
    def addMeson(self,meson):
        """Add meson spectroscopy."""
        self.cycle[self.cycnt].meson.append(meson)
        return meson
    def addBaryon(self,baryon):
        """Add baryon spectroscopy (*NOTE* unimplemented)."""
        raise
    def _bind_indices(self):
        """Bind indices to sources, propagators and quarks."""
        for cy in self.cycle:
            nbs = len(cy.bsource)
            for idx, s in zip(range(nbs),cy.bsource):
                s.id = idx
                pass
            nms = len(cy.msource)
            for idx, s in zip(range(nbs,nbs+nms),cy.msource):
                s.id = idx
                pass
            idx = 0
            for ps in cy.pset:
                for p in ps:
                    p.id = idx
                    idx += 1
                    pass
                pass
            for idx, q in zip(range(len(cy.quark)),cy.quark):
                q.id = idx
                pass
        return
    def generate(self,ostream=sys.stdout):
        """Write MILC prompts to output stream 'ostream'."""
        self._bind_indices()
        self.geometry.generate(ostream)
        for x in self.cycle:
            x.generate(ostream)
            pass
        return
    def dataflow(self):
        self._bind_indices()
        w = list()
        for x in self.cycle:
            w.append(x.dataflow())
            pass
        dinfo = { 'title': self.pName, 'workflow': w }
        return dinfo
    pass

class _Cycle_ks_spectrum:
    def __init__(self):
        self.gauge = None
        self.eigen = None
        self.bsource = list()
        self.msource = list()
        self.nextpropidx = 0
        self.pset = list()
        self.quark = list()
        self.meson = list()
        self.baryon = list()
        return
    def generate(self,ostream):
        self.gauge.generate(ostream)
        self.eigen.generate(ostream)
        print>>ostream, '#== PBP Masses =='
        print>>ostream
        print>>ostream, 'number_of_pbp_masses 0'
        print>>ostream
        print>>ostream, '#== Base Sources =='
        print>>ostream
        print>>ostream, 'number_of_base_sources', len(self.bsource)
        for x in self.bsource:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '#== Modified Sources =='
        print>>ostream
        print>>ostream, 'number_of_modified_sources', len(self.msource)
        for x in self.msource:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '#== KSsolveSets =='
        print>>ostream
        print>>ostream, 'number_of_sets', len(self.pset)
        for x in self.pset:
            x.generate(ostream)
        print>>ostream
        print>>ostream, '#== Quarks =='
        print>>ostream
        print>>ostream, 'number_of_quarks', len(self.quark)
        for x in self.quark:
            x.generate(ostream)
            pass
        print>>ostream
        print>>ostream, '#== Mesons =='
        print>>ostream
        print>>ostream, 'number_of_mesons', len(self.meson)
        for x in self.meson:
            x.generate(ostream)
        print>>ostream
        print>>ostream, '#== Baryons =='
        print>>ostream
        print>>ostream, 'number_of_baryons', len(self.baryon)
        for x in self.baryon:
            x.generate(ostream)
        return
    def dataflow(self):
        dinfo = list()
        dinfo.append(self.gauge.dataflow())
        for x in self.bsource:
            dinfo.append(x.dataflow())
            pass
        for x in self.msource:
            dinfo.append(x.dataflow())
            pass
        for x in self.pset:
            for y in x.dataflow():
                dinfo.append(y)
                pass
            pass
        for x in self.quark:
            dinfo.append(x.dataflow())
            pass
        for x in self.meson:
            dinfo.append(x.dataflow())
            pass
        for x in self.baryon:
            dinfo.append(x.dataflow())
            pass
        return dinfo
    pass
