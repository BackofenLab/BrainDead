def eval_rna_structure(sequence, dotbracket):
    from subprocess import Popen, PIPE
    cmd = " ".join(['RNAeval '])
    # Verify there si no whitespace in the strings
    assert len(sequence.split()) == 1 and len(dotbracket.split()) == 1
    if len(sequence) != len(dotbracket):
        raise RuntimeError("Mismatch length of sequence and structures: {} {}".format(sequence, dotbracket))
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    input_str = "{}\n{}\n".format(sequence, dotbracket)
    out, err = p.communicate(input=input_str)
    has_error = False
    if err:
        print "Error in calling eval_rna_structure \n", out, err
        has_error = True

    splits = [l for l in out.split("\n") if len(l.strip()) != 0]  # Remove empty lines
    assert len(splits) == 2
    line = splits[1]
    import re
    rnaeval_re = re.compile(r'([\.\(\)]+)\s+\(\s*([\d\.\-]+)\s*\)')  # (r'(['+ re.escape('.()')+']+)\s+(\[\s+\])')
    matches = rnaeval_re.match(line)
    if matches is None:
        raise RuntimeError("Error cannot match re for:\n".format(line))
        has_error = True
    else:
        groups = matches.groups()
    assert len(groups) == 2
    struct = groups[0]
    assert struct == dotbracket
    energy = (float)(groups[1])
#     print({'energy':energy, 'struct':struct})
    return energy, has_error

def parse_subpot_out(f, maxnum_subopts=None):
    import re
    subopt_re = re.compile(r'([\.\(\)]+)\s+\[?([\d\.\-]+)\]?')  #(r'(['+ re.escape('.()')+']+)\s+(\[\s+\])')

    suboptimals = list()
    with open(f) as in_file:
        for l, line in enumerate(in_file):
            if maxnum_subopts is not None and l > maxnum_subopts+1:
                print "Stopped reading subopts exceeding the maxnum"
                break
            if line.startswith(">") or re.match('[ACGUNacgu]+', line):
                continue
            line = line.strip()
            #print line
            matches = subopt_re.match(line).groups()
            assert len(matches) == 2
            struct = matches[0]
            energy = (float)(matches[1])
            suboptimals.append({'energy': energy, 'struct': struct})
    print "#suboptimal:", len(suboptimals)
    return suboptimals


def call_RNAsubopt(sequence, seq_name, deltaEnergy=3, out_dir='./', use_cache=False):
    '''Runs Vienna RNAsubopt for the given sequence
    '''
    import os
    outfile = out_dir+'/{}.subopt.out'.format(seq_name)
    if use_cache is True and os.path.isfile(outfile):
        print "call_RNAsubopt re-using cache"
        parse_subpot_out
        return outfile

    from subprocess import Popen, PIPE

    cmd = " ".join(['RNAsubopt', '-e {}'.format(deltaEnergy), '-s '])
    # print cmd
    # Verify there si no whitespace in the strings
    assert len(sequence.split()) == 1

    out_handle = open(outfile, 'w')
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=out_handle, stderr=PIPE)
    input_str = ">{}\n{}\n".format(seq_name, sequence)
    # print input_str
    out, err = p.communicate(input=input_str)
    if err:
        raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
    out_handle.close()

    return outfile


def call_RNAfold(sequence, seq_name):
    '''Runs Vienna RNAfold for MFE prediction, for all sequences inside input fasta file
    '''
    from subprocess import Popen, PIPE
    cmd = " ".join(['RNAfold '])
    # Verify there si no whitespace in the strings
    assert len(sequence.split()) == 1
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    input_str = ">{}\n{}\n".format(seq_name, sequence)
    out, err = p.communicate(input=input_str)
    has_error = False
    if err:
        raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
#     print out
    lines = [l for l in out.split("\n") if len(l.strip()) != 0]  # Remove empty lines
    print lines
    assert len(lines) == 3

    line = lines[2]
    import re
    rnaeval_re = re.compile(r'([\.\(\)]+)\s+\(\s*([\d\.\-]+)\s*\)')  # (r'(['+ re.escape('.()')+']+)\s+(\[\s+\])')
    matches = rnaeval_re.match(line)
    if matches is None:
        raise RuntimeError("Error cannot match re for:\n".format(line))
        has_error = True
    else:
        groups = matches.groups()
    assert len(groups) == 2
    struct = groups[0]
    assert len(struct) == len(sequence)
    mfe_energy = (float)(groups[1])
#     print({'energy':energy, 'struct':struct})
    return struct, mfe_energy, has_error


def call_RNAfold_pf(in_seq, seq_name):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file
    Returns:  [mfe_struct, mfe_energy, pf_energy] '''

    from subprocess import Popen, PIPE
    RNAFOLD = 'RNAfold -p2 '
    p = Popen(('echo "%s" | ' % in_seq) + RNAFOLD, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)

    out, err = p.communicate()
    if err:
        print "Error in calling RNAfold for ", in_seq
        print out
        print err
        raise RuntimeError
    lines = out.split('\n')
    #     print out
    mfe_line = lines[1]  # mfe line
    assert len(mfe_line.split()) >= 2
    mfe_struct = mfe_line.split()[0]
    assert len(mfe_struct) == len(in_seq)
    mfe_energy = float(mfe_line.split()[-1].replace(')', '').replace('(', ''))

    pf_line = lines[2]  # pf line
    assert len(pf_line.split()) >= 2
    pf_energy = float(pf_line.split()[-1].replace(']', '').replace('[', ''))

    import re
    my_re = re.compile(r'\s*frequency of mfe structure in ensemble\s+([\d\.\-e]+);\sensemble diversity\s+([\d\.\-]+)')
    freq_line = lines[4]
    match = my_re.match(freq_line)
    if match is None or match.groups() is None:
        print "Error unexpected frequency line format, found:\n  {}\n  {}\n".format(in_seq, freq_line)
        raise RuntimeError
    mfe_prob = (float)(match.groups()[0])
    diversity = (float)(match.groups()[1])
    return mfe_struct, mfe_energy, pf_energy  # , mfe_prob, diversity

