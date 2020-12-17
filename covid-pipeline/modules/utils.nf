process concat_elements_to_single_string{
    input:
      val string_value_list

    output:
      val concatenated_string

    script:
      concatenated_string = string_value_list.join(" ")

    """
    """
}