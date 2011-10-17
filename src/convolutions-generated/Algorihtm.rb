module ConvolutionGenerator
  class Variable
    attr_reader :name, :type, :policy
    def initialize(name, type, policy)
      @name = name
      @type = type
      @policy = policy
    end
    def to_s
       return @name
    end
  end

  class Algorithm
    attr_reader :name, :variables, :parameters
    def initialize(name, variables, parameters)
      @name = name
      @variables = variables
      @parameters = parameters
    end  
  end


  class For < Algorithm
    attr_reader :start, :stop, :increment, :body
    def initialize(variable, start, stop, increment, body)
      super("For", [ variable ], nil)
      @start = start
      @stop = stop
      @body = body
      @increment = increment
    end

    def to_c
      res = "for( #{@variables.first.to_s} = #{@start.to_s} ; #{@variables.first.to_s} <= #{@stop.to_s} ; #{@variables.first.to_s} += #{@increment} ) {\n"
      res += body.to_s
      res += "}\n"
      return res
    end

    def to_fortran
      res = "do #{@variables.first.to_s} = #{@start.to_s}, #{@stop.to_s}, #{@increment}\n"
      res += body.to_s
      res += "enddo\n"
      return res
    end

  end
end
