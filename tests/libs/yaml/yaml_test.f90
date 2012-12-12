program toto
   use yaml_output
   call yaml_open_map("Essai")
      call yaml_map("J'ai une grande phrase qui ne veut rien dire mais qui dépasse la taille",.true.)
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("un",1)
      call yaml_map("deux",2)
      call yaml_close_map()
      call yaml_open_map("toto",flow=.true.)
      call yaml_map("un",1)
      call yaml_map("deux",2)
      call yaml_close_map()
   call yaml_close_map()
end program toto
