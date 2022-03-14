from collections import OrderedDict
from sympy import Symbol, Function

def collect_recursively(expr, syms, func=None, evaluate=True, exact=False, distribute_order_term=True):
    """
    Use `sympy.collect` function recursively.

    This function behaves just like the standard sympy `collect`
    function, except that `syms` is now a list of tuples.  This
    function steps through that list, performing `collect` with each
    tuple, and then performing `collect` with successive tuples on
    each of the collected terms.

    The optional arguments are the same as the usual function,
    except that func is only applied at the innermost level.

    """
    from sympy import collect
    if len(syms)==0:
        return collect(expr) # should raise an exception
    def collect_factory(syms, func, evaluate, exact, distribute_order_term):
        return lambda x: collect(x, syms, func, evaluate, exact, distribute_order_term)
    syms[-1] = collect_factory(syms[-1], func, evaluate, exact, distribute_order_term)
    for i in range(len(syms)-2, -1, -1):
        syms[i] = collect_factory(syms[i], syms[i+1], evaluate, exact, distribute_order_term)
    return syms[0](expr)

class PNSymbol(Symbol) :
    """This is the basic object created by calls to `AddVariable`, etc.,
    and is a simple subclass of sympy.Symbol.

    The method `__new__` is always called first, since this is an
    immutable object, which creates the object, allocating memory for
    it.  Since `__new__` actually returns an object, `__init__` is
    then called with the same arguments as `__new__`.  This is why we
    throw away the custom arguments in `__new__`, and throw away the
    rest in `__init__`.

    """
    def __new__(cls, name, constant=None, fundamental=None, substitution=None, substitution_atoms=None, datatype=None, **assumptions) :
        from sympy import Symbol
        return Symbol.__new__(cls, name, **assumptions)
    def __init__(self, name, constant=None, fundamental=None, substitution=None, substitution_atoms=None, datatype=None, **kwargs) :
        if not fundamental and isinstance(substitution, str) and not substitution_atoms:
            raise ValueError('Either `substitution` must be a sympy expression, '
                             +'or `substitution_atoms` must be non-empty for derived quantities.')
        self.constant = constant
        self.fundamental = fundamental
        self.substitution = substitution
        if substitution_atoms:
            self.substitution_atoms = substitution_atoms
        else:
            try:
                self.substitution_atoms = self.substitution.atoms(Symbol)
            except AttributeError:
                self.substitution_atoms = None
        self.datatype = datatype
    def ccode(self, **args):
        from sympy import ccode, N#, horner
        if self.fundamental:
            return str(self)
        if hasattr(self.substitution, '__iter__'): # Check only for __iter__ so that strings don't get caught
            return '{' + ', '.join([ccode(x, **args) for x in self.substitution]) + '}'
        if isinstance(self.substitution, str):
            return self.substitution
        code = self.substitution
        try:
            code=N(code)
        except:
            pass
        # try:
        #     code=horner(code, wrt=args.pop('wrt', None))
        # except:
        #     pass
        return ccode(code, **args)
    def pycode(self, **args):
        from sympy import pycode, N#, horner
        if self.fundamental:
            return str(self)
        if hasattr(self.substitution, '__iter__'): # Check only for __iter__ so that strings don't get caught
            return '{' + ', '.join([pycode(x, **args) for x in self.substitution]) + '}'
        if isinstance(self.substitution, str):
            return self.substitution
        code = self.substitution
        try:
            code=N(code)
        except:
            pass
        # try:
        #     code=horner(code, wrt=args.pop('wrt', None))
        # except:
        #     pass
        return pycode(code, **args)

class PNCollection(OrderedDict) : # subclass of OrderedDict
    """Subclass of `OrderedDict` to hold PN variables, each of which is a
    subclass sympy `Symbol`

    This class has a few extra functions to add variables nicely, and
    include them in the calling scope.  The PN variables are just like
    sympy `Symbol`s, except they also remember these things:

      * `constant` (boolean: doesn't need to be updated once it's defined)
      * `fundamental` (boolean: not defined in terms of other things)
      * `substitution` (string or sympy expression: re-express non-fundamental object in terms of fundamentals)
      * `substitution_atoms` (atomic sympy objects in terms of which the variable is defined)
      * `datatype` (optional name for datatype of variable in code output;
                    if None, the basic real datatype of the output language is assumed)

    """
    def __init__(self, *args):
        OrderedDict.__init__(self, *args)
    def _AddVariable(self, name, **args) :
        from inspect import currentframe
        from sympy import Basic, FunctionClass
        add_to_globals = args.pop('add_to_globals', True)
        frame = currentframe().f_back.f_back
        if not args['fundamental'] and isinstance(args['substitution'], str) and not args['substitution_atoms']:
            # Try to add appropriate PNSymbols to list of atoms in substitution
            from sympy import sympify
            args['substitution_atoms'] = list(sympify(args['substitution']).atoms(Symbol))
            atom_strings = [str(a) for a in args['substitution_atoms']]
            for k,v in self.items():
                if v in atom_strings:
                    args['substitution_atoms'][atom_strings.index(v)] = k
        try:
            sym = PNSymbol(name, **args)
            if sym is not None and add_to_globals:
                if isinstance(sym, Basic):
                    frame.f_globals[sym.name] = sym
                elif isinstance(sym, FunctionClass):
                    frame.f_globals[sym.__name__] = sym
        finally:
            del frame
        self.update({sym:name})
        return sym
    def AddVariable(self, name, constant=False, fundamental=False, substitution=None, substitution_atoms=None, datatype=None, **args) :
        return self._AddVariable(name, constant=constant, fundamental=fundamental,
                                 substitution=substitution, substitution_atoms=substitution_atoms, datatype=datatype, **args)
    def AddBasicConstants(self, names, datatype=None, **args) :
        from re import split
        names = split(',| ', names)
        args['constant'] = True
        args['fundamental'] = True
        args['substitution'] = args.pop('substitution', None)
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = datatype
        for name in names :
            if name :
                self._AddVariable(name, **args)
    def AddBasicVariables(self, names, **args) :
        from re import split
        names = split(',| ', names)
        args['constant'] = False
        args['fundamental'] = True
        args['substitution'] = args.pop('substitution', None)
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        for name in names :
            if name :
                self._AddVariable(name, **args)
    def AddDerivedConstant(self, name, substitution, **args) :
        args['constant'] = True
        args['fundamental'] = False
        args['substitution'] = substitution
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        self._AddVariable(name, **args)
    def AddDerivedVariable(self, name, substitution, **args) :
        args['constant'] = False
        args['fundamental'] = False
        args['substitution'] = substitution
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        self._AddVariable(name, **args)
PNCollection.AddExpression = PNCollection.AddDerivedVariable
