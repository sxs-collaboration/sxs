class CodeConstructor:
    """Contains lists of variables and expressions to be written as code.

    `CodeConstructor` objects contain:
      1) An ordered list of atoms for the code to use
      2) A PNCollection of PNSymbol objects
      3) A PNCollection of expressions to be calculated

    Once the `CodeConstructor` is initialized with these objects, it
    can be used to construct various types of code.  For example, the
    `CppDeclarations` method will output a list of declarations of the
    atoms.  Similar methods are available for function input
    arguments, class initializer lists, and the final evaluations
    needed to calculate the input `Expressions`.

    Support for other languages or constructions can be added by
    adding more method functions to this class.

    Note that it is generally necessary to obey a strict ordering for
    defining variables.  The functions of this class assume that the
    ordering in which the variables were defined in python should
    remain the same in the output code.  Because of the structure of
    the `PNCollection` objects, it should be hard to define the
    variables out of order, so this should not require anything from
    the user.  (That is why `PNCollection` is a subclass of the basic
    `OrderedDictionary` object.)  However, if new functions are added
    here, they must obey that ordering.

    """
    def __init__(self, Variables, Expressions,NumTerms=None,DenTerms=None):
        AtomSet = set([])
        self.Variables = Variables
        self.Expressions = Expressions
        for Expression in self.Expressions:
            try:
                AtomSet.update(Expression.substitution_atoms)
            except TypeError:
                pass
        LastAtomsLength = 0
        while(len(AtomSet) != LastAtomsLength):
            LastAtomsLength = len(AtomSet)
            for Atom in list(AtomSet): 
                try:
                    if (Atom.substitution_atoms):
                        AtomSet.update(Atom.substitution_atoms)
                except:
                    try:
                        Atom=NumTerms.get(str(Atom))
                        AtomSet.update(Atom.substitution_atoms)
                    except:
                        try:
                            Atom=DenTerms.get(str(Atom))
                            AtomSet.update(Atom.substitution_atoms)
                        except:
                            pass
        self.Atoms = []
        for sym in self.Variables:
            if sym in AtomSet or str(sym) in set(['xHat','yHat','zHat','S_chi1','S_chi2','R_S1','R_S2','rfrak_chi1','rfrak_chi2','rfrak_frame']):
                self.Atoms.append(sym)

    @staticmethod
    def const(e):
        if e.constant:
            return 'const '
        else:
            return ''

    @staticmethod
    def dtype(e):
        if e.datatype:
            return e.datatype
        else:
            return 'float64'

    def AddDependencies(self, Expressions):
        AtomSet = set([])
        for Expression in Expressions:
            if (Expression.substitution_atoms):
                AtomSet.update(Expression.substitution_atoms)
        LastAtomsLength = 0
        while(len(AtomSet) != LastAtomsLength):
            LastAtomsLength = len(AtomSet)
            for Atom in list(AtomSet):
                if (Atom.substitution_atoms):
                    AtomSet.update(Atom.substitution_atoms)
        OldAtoms = self.Atoms[:]
        self.Atoms = []
        for sym in self.Variables:
            if sym in AtomSet or sym in OldAtoms:
                self.Atoms.append(sym)
        # for Expression in Expressions:
        #     try:
        #         for sym in Expression.Variables:
        #             if sym in AtomSet or sym in OldAtoms:
        #                 self.Atoms.append(sym)
        #     except:
        #         pass

    def CppDeclarations(self, Indent=4):
        """Create declaration statements for C++

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the declaration
        list should be

            const double m1, m2;
            double t, x;

        The code knows which atoms need to be declared at the
        beginning, and which ones should be `const`, for example.  For
        C++, the default datatype is `double`; if the atom was created
        with a different datatype, that will be used appropriately.

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ''
        wrapper.subsequent_indent = ''
        datatype = ''
        Declarations = ''
        Names = []
        for atom in self.Atoms:
            thisdatatype = CodeConstructor.const(atom) + CodeConstructor.dtype(atom) + ' '
            if thisdatatype != datatype:
                if Names:
                    Declarations += wrapper.fill(', '.join(Names)) + ";\n"
                Names = []
                datatype = thisdatatype
                wrapper.initial_indent = ' '*Indent + thisdatatype
                wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
            Names.append(self.Variables[atom])
        if Names:
            Declarations += wrapper.fill(', '.join(Names)) + ";\n"
        return Declarations.rstrip()

    def CppInputArguments(self, Indent=12):
        """Create basic input arguments for C++

        The fundamental variables are listed, along with their data
        types and `const` if the variable is constant.  This would be
        an appropriate string to represent the input arguments for a
        function or class constructor to calculate the `Expressions`
        of this CodeConstructor object.

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the input argument
        list should be

            const double m1_i, const double m2_i, double t_i, double x_i

        """
        #from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent
        InputArguments = ['{0}_i'.format(self.Variables[atom])
                          for atom in self.Atoms if atom.fundamental and str(atom).find('rfrak_chi')==-1]
        return ', '.join(InputArguments)

         
    
    def CppInitializations(self, Indent=4):
        """Create initialization list for C++

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the initialization
        list should be

            m1(m1_i), m2(m2_i), t(t_i), x(x_i)

        The quantities m1_i, etc., appear in the input-argument list
        output by the method `CppInputArguments`.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent
        def Initialize(atom):
            A=str(atom.substitution)
            if A[0:3]=='dot':
                A='np.array(['+A+'])'
            k=0
            for i in self.Atoms:
                k+=A.find(str(i))
                if i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
            if k==-len(self.Atoms):
                A='np.array(['+A+'])'
            return A
        def Initialization(atom):
            if atom.constant:
                if atom.datatype and (atom.datatype=='double'):
                    return '    Cons.{0}={1}'.format(self.Variables[atom], len(atom.substitution))
                if atom.fundamental and str(self.Variables[atom]).find('Hat')==-1 and str(self.Variables[atom]).find('S')==-1:
                    return '    Cons.{0}=np.array([{0}_i])'.format(self.Variables[atom])
                elif atom.fundamental and (str(self.Variables[atom]).find('Hat')!=-1 or str(self.Variables[atom]).find('S')!=-1):
                    return '    Cons.{0}={0}_i'.format(self.Variables[atom])
                else:
                    return '    Cons.{0}={1}'.format(self.Variables[atom], Initialize(atom))
            else:
                if atom.datatype and (atom.datatype=='double'):
                    return '   {0}={1}'.format(self.Variables[atom], len(atom.substitution))
                if atom.fundamental and str(atom).find('rfrak_chi')==-1:
                    return '    {0}=np.array([{0}_i])'.format(self.Variables[atom])
                elif atom.fundamental and str(atom).find('rfrak_chi')!=-1:
                    return '    {0}=np.array([0.0,0.0])'.format(self.Variables[atom])
                else:
                    return '    {0}={1}'.format(self.Variables[atom], Initialize(atom))
        Initializations  = [Initialization(atom) for atom in self.Atoms if atom.constant or str(atom).find('chiVec')!=-1 or str(atom).find('R_S')!=-1 or str(atom).find('rfrak_chi')!=-1]
        return '\n'.join(Initializations)
    
    def CppInitializationModes(self, Indent=4):
        """Create initialization list for C++

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the initialization
        list should be

            m1(m1_i), m2(m2_i), t(t_i), x(x_i)

        The quantities m1_i, etc., appear in the input-argument list
        output by the method `CppInputArguments`.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent
        def Initialize(atom):
            A=str(atom.substitution)
            if A[0:3]=='dot':
                A='np.array(['+A+'])'
            if str(atom).find('hHat')!=-1 and A.find('I')==-1:
                A=A+' + 0.0*I'
            k=0
            for i in self.Atoms:
                k+=A.find(str(i))
                if i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
            if k==-len(self.Atoms):
                A='np.array(['+A+'])'
            return A
        def Initialization(atom):
            if atom.constant:
                if atom.datatype and (atom.datatype=='double'):
                    return '    Cons.{0}={1}'.format(self.Variables[atom], len(atom.substitution))
                if atom.fundamental and str(self.Variables[atom]).find('Hat')==-1 and str(self.Variables[atom]).find('S')==-1:
                    return '    Cons.{0}=np.array([{0}_i])'.format(self.Variables[atom])
                elif atom.fundamental and (str(self.Variables[atom]).find('Hat')!=-1 or str(self.Variables[atom]).find('S')!=-1):
                    return '    Cons.{0}={0}_i'.format(self.Variables[atom])
                else:
                    return '    Cons.{0}={1}'.format(self.Variables[atom], Initialize(atom))
        Initializations  = [Initialization(atom) for atom in self.Atoms if atom.constant]
        return '\n'.join(Initializations)
    
    def CppRead(self):
        return '\n'.join(['        {0}=self.{0}'.format(self.Variables[atom]) for atom in self.Atoms])

    def ConsArgu(self):
        return ','.join(['{0}'.format(self.Variables[atom]) for atom in self.Atoms if atom.constant])
    
    def ConsSpec(self):
        def spec(atom):
            A=str(atom)
            if A.find('hHat')==-1:
                return '(\'{0}\', float64[:])'.format(self.Variables[atom])
            else:
                return '(\'{0}\', complex128[:])'.format(self.Variables[atom])
        return ','.join([spec(atom) for atom in self.Atoms if atom.constant])
    
    def ConsInit(self):
        return '\n'.join(['        self.{0}={0}'.format(self.Variables[atom]) for atom in self.Atoms if atom.constant])
    
    def VarsArgu(self):
        return ','.join(['{0}'.format(self.Variables[atom]) for atom in self.Atoms if not atom.constant])
    
    def VarsSpec(self):
        def spec(atom):
            A=str(atom)
            if A.find('hHat')==-1:
                return '(\'{0}\', float64[:])'.format(self.Variables[atom])
            else:
                return '(\'{0}\', complex128[:])'.format(self.Variables[atom])
        return ','.join([spec(atom) for atom in self.Atoms if not atom.constant])
    
    def VarsInit(self):
        return '\n'.join(['        self.{0}={0}'.format(self.Variables[atom]) for atom in self.Atoms if not atom.constant])
    
    def CppInputSpecs(self):
        return ','.join(['("{0}_i", float64)'.format(self.Variables[atom]) for atom in self.Atoms if atom.fundamental])
        

    def CppEvaluations(self, Indent=4, GlobalDeclaration=False):
        """Evaluate all derived variables in C++

        This function uses the `substitution` expressions for the
        derived variables.  This output is appropriate for updating
        the values of the variables at each step of an integration,
        for example.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent + '  '
        def Evaluate(atom):
            A=str(atom.substitution)
            if A[0:3]=='dot':
                A='np.array(['+A+'])'
            for i in self.Atoms:
                if i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
                elif not i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
            return A
        
        def Evaluation(atom):
            def Ccode(a) :
                try:
                    return a.pycode()
                except :
                    from sympy.printing import pycode
                    return pycode(a)
            if atom.datatype and atom.datatype=='double':
                return '\n'.join(['    {0}[{1}] = {2};'.format(self.Variables[atom], i, Ccode(atom.substitution[i]))
                                  for i in range(len(atom.substitution))])
            else:
                return '    Vars.{0} = {1}'.format(self.Variables[atom], Evaluate(atom))
        ReturnExpression='\n'.join([Evaluation(atom) for atom in self.Atoms if not atom.fundamental and not atom.constant])
        return ReturnExpression
    
    def CppEvaluationModes(self, Indent=4, GlobalDeclaration=False):
        """Evaluate all derived variables in C++

        This function uses the `substitution` expressions for the
        derived variables.  This output is appropriate for updating
        the values of the variables at each step of an integration,
        for example.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent + '  '
        def Evaluate(atom):
            A=str(atom.substitution)
            if A[0:3]=='dot':
                A='np.array(['+A+'])'
            if str(atom).find('hHat')!=-1 and A.find('I')==-1:
                A=A+' + 0.0*I'
            for i in self.Atoms:
                if i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
                elif not i.constant:
                    j=0
                    while j!=-1 and j+len(str(i))<=len(A):
                        j=A.find(str(i),j)
                        if j!=-1:
                            if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                A=B
                                j=j+5+len(str(i))
                            else:
                                j=j+1
            return A
        
        def Evaluation(atom):
            def Ccode(a) :
                try:
                    return a.pycode()
                except :
                    from sympy.printing import pycode
                    return pycode(a)
            if atom.datatype and atom.datatype=='double':
                return '\n'.join(['    {0}[{1}] = {2};'.format(self.Variables[atom], i, Ccode(atom.substitution[i]))
                                  for i in range(len(atom.substitution))])
            else:
                return '    Vars.{0} = {1}'.format(self.Variables[atom], Evaluate(atom))
        ReturnExpression='\n'.join([Evaluation(atom) for atom in self.Atoms if not atom.fundamental and not atom.constant])
        return ReturnExpression

    def CppEvaluateExpressions(self, Indent=4, Expressions=None):
        """Declare and define the `Expressions` for C++

        The output of this function declares are defines the
        `Expressions` as individual variables.  An optional dictionary
        of expressions allows just a subset of this object's
        expressions to be output; if this argument is not present, all
        will be output.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent+'  '
        def Express(atom):
            A=str(atom)
            for i in self.Atoms:
                try:
                    if i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                    A=B
                                    j=j+5+len(str(i))
                                else:
                                    j=j+1
                    elif not i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    if j==0:
                                        B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                        A=B
                                        j=j+5+len(str(i))
                                    elif not (A[j-1].isalnum() or A[j-1]=='_'):
                                        B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                        A=B
                                        j=j+5+len(str(i))
                                    else:
                                        j=j+1
                                else:
                                    j=j+1
                except:
                    pass
            return A
        
        Evaluations = []
        if not Expressions:
            Expressions=self.Expressions
        for Expression in Expressions:
            try:
                Evaluations.append('    {0} = {1}'.format(Expressions[Expression], Express(Expression.pycode())))
            except TypeError:
                pass
        return '\n'.join(Evaluations)
    
    def CppEvaluateT45(self, NumTerms, DenTerms, Constructor, Expressions=None):
        """Declare and define the `Expressions` for C++

        The output of this function declares are defines the
        `Expressions` as individual variables.  An optional dictionary
        of expressions allows just a subset of this object's
        expressions to be output; if this argument is not present, all
        will be output.

        """
        from textwrap import TextWrapper
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent
        #wrapper.subsequent_indent = wrapper.initial_indent+'  '
        def Express(atom):
            A=str(atom)
            while (A.find('Num')!=-1):
                j=A.find('Num')
                Subs=str(NumTerms.get(str(A[j:j+4])))
                k=Subs.find('0*')
                B=''
                if k!=0 and(A[k-1].isalnum() or A[k-1]=='_' or A[k-1]!='.' or k==-1):
                    B=A[0:j]+Subs+A[j+4:]
                else:
                    B=A[0:j]+'0'+A[j+4:]
                A=B
            while (A.find('Den')!=-1):
                j=A.find('Den')
                Subs=str(DenTerms.get(str(A[j:j+4])))
                k=Subs.find('0*')
                B=''
                if k!=0 and(A[k-1].isalnum() or A[k-1]=='_' or A[k-1]!='.' or k==-1):
                    B=A[0:j]+Subs+A[j+4:]
                else:
                    B=A[0:j]+'0'+A[j+4:]
                A=B
            while (A.find('PolynomialVariable')!=-1):
                j=A.find('PolynomialVariable')
                B=A[0:j]+'v'+A[j+18:]
                A=B
            for i in Constructor.Atoms:
                try:
                    if i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                    A=B
                                    j=j+5+len(str(i))
                                else:
                                    j=j+1
                    elif not i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    if j==0:
                                        B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                        A=B
                                        j=j+5+len(str(i))
                                    elif not (A[j-1].isalnum() or A[j-1]=='_'):
                                        B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                        A=B
                                        j=j+5+len(str(i))
                                    else:
                                        j=j+1
                                else:
                                    j=j+1
                except:
                    pass
            return A
        
        Evaluations = []
        if not Expressions:
            Expressions=self.Expressions
        for Expression in Expressions:
            try:
                Evaluations.append('    {0} = {1}'.format(Expressions[Expression], Express(Expression.pycode())))
            except TypeError:
                pass
        return '\n'.join(Evaluations)

    def CppExpressionsAsFunctions(self, Indent=4, OrderString='', Expressions=None):
        """Define functions to calculate the `Expressions` in C++

        The output of this function gives C++ functions to calculate
        the `Expressions`, assuming the functions are member methods
        in a class, and so can access the atoms of the expression
        without explicit arguments.  An optional dictionary of
        expressions allows just a subset of this object's expressions
        to be output; if this argument is not present, all will be
        output.

        """
        def Express(atom):
            A=str(atom)
            for i in self.Atoms:
                try:
                    if i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    B=A[0:j]+'Cons.{0}'.format(i)+A[j+len(str(i)):]
                                    A=B
                                    j=j+5+len(str(i))
                                else:
                                    j=j+1
                    elif not i.constant:
                        j=0
                        while j!=-1 and j+len(str(i))<=len(A):
                            j=A.find(str(i),j)
                            if j!=-1:
                                if j+len(str(i)) ==len(A) or not (A[j+len(str(i))].isalnum() or A[j+len(str(i))]=='_'):
                                    B=A[0:j]+'Vars.{0}'.format(i)+A[j+len(str(i)):]
                                    A=B
                                    j=j+5+len(str(i))
                                else:
                                    j=j+1
                except:
                    pass
            return A
        
        def dtype(e):
            if e.datatype:
                return e.datatype
            else:
                return 'double'
        from textwrap import TextWrapper
        from Utilities.PNObjects import PNCollection
        #wrapper = TextWrapper(width=120)
        #wrapper.initial_indent = ' '*Indent + '  return '
        #wrapper.subsequent_indent = ' '*Indent + '    '
        Evaluations = []
        if not Expressions:
            Expressions=self.Expressions
        for Expression in Expressions:
            ExprColl = PNCollection()
            for atom in Expression.substitution_atoms:
                if atom not in self.Variables:
                    try:
                        ExprColl.AddDerivedVariable(str(atom), atom.substitution,
                                                    substitution_atoms=atom.substitution_atoms,
                                                    datatype=atom.datatype)
                    except TypeError:
                        pass
            MiniConstructor = CodeConstructor(self.Variables, ExprColl)
            if str(Expression) != 'OrbitalAngularMomentum':
                Evaluations.append(
                    '@njit\n'+'def ' + Expressions[Expression] + '_{0}(Cons,Vars):\n'.format(OrderString)
                    + MiniConstructor.CppEvaluateExpressions(Indent+2) + '\n    return '
                    + Express(Expression.pycode()) + '\n'
            )
        return '\n'.join(Evaluations)