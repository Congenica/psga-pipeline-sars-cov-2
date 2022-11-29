/*
 * Main workflow for the pathogen: Staphylococcus aureus.
 */
params.str = 'Hello world!'

process splitLetters {
  output:
    path 'chunk_*'

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}


workflow psga {
  splitLetters | flatten | convertToUpper | view { it.trim() }
}
