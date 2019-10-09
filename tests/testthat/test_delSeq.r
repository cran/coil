
test_that("A sequence with a deletion is framed and translated properly", {

  sequence = 'accctttacttaatctttggtgcatgagcaggaatagtaggtacagcccttagcttgcttattcgagcagaattaagccaacctggcacactcctgggagacgatcagatctacaatgttatcgtaactgctcacgcttttgtaataattttttttatggttatacctgtaataattggtgggttcggaaactgattagtgcctttaataattggtgcaccggacatagctttcccacgaataaataacataagcttttgactgctacccccctccctcctattacttttggcctctgctggagttgaagccggagccggaactggttgaacagtttatccccccctcgcaagtaatatagcccacgctggggcatcagtagacttagctattttctcgctccatttagcggtatttcctcaattcttgcctctatcaactttattacaaccattattaatataaaaccgcctgccatctctcaatatcaaacacccctttttgtttgatctattcttgtaaccacagtcctactcctcctttcacttcctgttcttgcagccgcaattacaatactacttaccgaccgtaatttaaacacaacattttttgatcctgctggtgggggtgacccaattctttaccaacattta'
  sequence_framed = 'accctttacttaatctttggtgcatgagcaggaatagtaggtacagcccttagcttgcttattcgagcagaattaagccaacctggcacactcctgggagacgatcagatctacaatgttatcgtaactgctcacgcttttgtaataattttttttatggttatacctgtaataattggtgggttcggaaactgattagtgcctttaataattggtgcaccggacatagctttcccacgaataaataacataagcttttgactgctacccccctccctcctattacttttggcctctgctggagttgaagccggagccggaactggttgaacagtttatccccccctcgcaagtaatatagcccacgctggggcatcagtagacttagctattttctcgctccatttagcggtatttcctcaattcttgcctctatcaactttattacaaccattattaatataaaaccgcctgccatctctcaatatcaaacacccctttttgtttgatctattcttgtaaccacagtcctactcctcctttcacttcctgttcttgcagccgcaattacaatactacttaccgaccgtaatttaaacacaacattttttgatcctgctggtgggggtgacccaattctttaccaacattta'

  sequence_AAcensored = "TLYLIFGAWAG?VGTALSLLIRAELSQPGTLLGDDQIYNVIVTAHAFV?IFFMV?PV?IGGFGNWLVPL?IGAPD?AFPR?NN?SFWLLPPSLLLLLASAGVEAGAGTGWTVYPPLASN?AHAGASVDLAIFSLHLAVFPQFLPLSTLLQPLL??NRLPSLNI?HPFLFDLFL?PQSYSSFHFLFLQPQLQYYLPTVI?TQHFLILLVGVTQFFTNI-"
  sequence_AA5 = "TLYLIFGAWAGMVGTALSLLIRAELSQPGTLLGDDQIYNVIVTAHAFVMIFFMVMPVMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSLLLLLASAGVEAGAGTGWTVYPPLASNMAHAGASVDLAIFSLHLAVFPQFLPLSTLLQPLLM*NRLPSLNIKHPFLFDLFL*PQSYSSFHFLFLQPQLQYYLPTVI*TQHFLILLVGVTQFFTNI"

  dat = coi5p(sequence )
  expect_equal(dat$raw, sequence)
  expect_identical(dat$name, character(0))

  dat = frame(dat)
  expect_equal(dat$framed, sequence_framed)

  dat = translate(dat)
  expect_equal(dat$aaSeq, sequence_AAcensored)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, TRUE)
  expect_equal(dat$stop_codons, FALSE)

  dat = translate(dat, trans_table = 5)
  expect_equal(dat$aaSeq, sequence_AA5)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, TRUE)
  expect_equal(dat$stop_codons, TRUE)

})
