class Dimension
  def initialize(name, multiples, notMultiple, paded, leadingDimension)
    @name = name
    @multiples = multiples
    @notMultiple = notMultiple
    @paded = paded
    @leadingDimension = leadingDimension
  end
end

class Data
  def initialize(typeSize, dimensions)
    @typeSize = typeSize
    @dimensions = dimensions
  end

end

class Pattern
  def initialize(name, inDimensions, outDimensions, vectorSize)
    @name = name
    @inDimensions = inDimensions
    @outDimensions = outDimenisons
    @vectorSize = vectorSize
  end
  
  def dimensionNumber
    return @outDimensions.length
  end

  def apply(dataIn, dataOut, baseResName, baseDataName, baseFilterName)
    raise "Incompatible dimensions!" if ( dataIn.dimensions.length =! @inDimensions.length or dataOut.dimensions.length =! @outDimensions.length )
    
  end

end
