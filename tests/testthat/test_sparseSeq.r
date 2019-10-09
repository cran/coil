
test_that("A sequence with missing information is framed and translated properly", {

  sequence = 'actttatatctaatttttggggcatgatcagcgatagttggaacagccttnnnnnnnnnnattcgagcagaactaggtcaacctggaagtttaattggagacgaccagatttataatgtgattgttactgcccacgcttttattataattttttttatggttatgcctattataattgggggattcggaaattgacttctacccttaatactaggagctccagatatagcttttccacgtctcaataatataagattttggttattaccccctgccttaatgttattaattagaggatccctagtagaagccggggcgggaactggatggacagtatacccacctctagccagcaatatcgcacnnnnnnnnnnctctgtagacctctctattttttcattacacttagcgggagcttcctctcttttaggagctatcaacttcatatccacagtaattaatatacgagctgaaactttaacattcgaccgtttacctttatttgtatgtagagtatttgtgacggtgattcttctattattatctttacctgtgttagccggagccattactatgcttttaactgaccggaatcttaatacnnnnnttttcgaccctacagggggaggagaccctattttataccaacactta'
  sequence_framed = 'actttatatctaatttttggggcatgatcagcgatagttggaacagccttnnnnnnnnnnattcgagcagaactaggtcaacctggaagtttaattggagacgaccagatttataatgtgattgttactgcccacgcttttattataattttttttatggttatgcctattataattgggggattcggaaattgacttctacccttaatactaggagctccagatatagcttttccacgtctcaataatataagattttggttattaccccctgccttaatgttattaattagaggatccctagtagaagccggggcgggaactggatggacagtatacccacctctagccagcaatatcgcacnnnnnnnnnnctctgtagacctctctattttttcattacacttagcgggagcttcctctcttttaggagctatcaacttcatatccacagtaattaatatacgagctgaaactttaacattcgaccgtttacctttatttgtatgtagagtatttgtgacggtgattcttctattattatctttacctgtgttagccggagccattactatgcttttaactgaccggaatcttaatacnnnnnttttcgaccctacagggggaggagaccctattttataccaacactta'
  sequence_AAcensored = "TLYLIFGAWSA?VGTA----IRAELGQPGSLIGDDQIYNVIVTAHAFI?IFFMVMPI?IGGFGNWLLPL?LGAPD?AFPRLNN??FWLLPPALMLLI?GSLVEAGAGTGWTVYPPLASNIA----SVDLSIFSLHLAGASSLLGAINF?STVIN?RAETLTFDRLPLFVC?VFVTVILLLLSLPVLAGAITMLLTDRNLN---FDPTGGGDPILYQHL"
  sequence_AA5 = "TLYLIFGAWSAMVGTA----IRAELGQPGSLIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLLPLMLGAPDMAFPRLNNMSFWLLPPALMLLISGSLVEAGAGTGWTVYPPLASNIA----SVDLSIFSLHLAGASSLLGAINFMSTVINMRAETLTFDRLPLFVCSVFVTVILLLLSLPVLAGAITMLLTDRNLNT--FDPTGGGDPILYQHL"

  dat = coi5p(sequence )
  expect_equal(dat$raw, sequence)
  expect_identical(dat$name, character(0))

  dat = frame(dat)
  expect_equal(dat$framed, sequence_framed)

  dat = translate(dat)
  expect_equal(dat$aaSeq, sequence_AAcensored)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

  dat = translate(dat, trans_table = 5)
  expect_equal(dat$aaSeq, sequence_AA5)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

})