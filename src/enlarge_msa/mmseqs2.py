from colabfold.colabfold import run_mmseqs2
import importlib_metadata
from pathlib import Path

AA_CLASSIC = "-XRHKDESTNQCGPAVILSMFYW"


def n_lower_chars(string):
    """Return the number of lowercase characters in a string.

    Parameters
    ----------
    string : str
        The string to count lowercase characters in.

    Returns
    -------
    int
        The number of lowercase characters in the string.
    """
    return sum(map(str.islower, string))


def only_insert(string):
    """Check if a string is only composed of '-' characters.

    Parameters
    ----------
    string : str
        The string to check.

    Returns
    -------
    bool
        True if the string is only composed of '-' characters, False otherwise.
    """

    for letter in string:
        if letter != "-":
            return False
    return True


def write_paired_alignement(f_out, a3m_lines_mmseqs2, start_indice, end_indice):
    """Write the paired alignement.

    Parameters
    ----------
    f_out : file
        The output file.
    a3m_lines_mmseqs2 : list of str
        The list of sequences.
    start_indice : list of int
        The list of start indices.
    end_indice : list of int
        The list of end indices.

    Returns
    -------
    None
    """

    # Paired alignement
    a3m_lines_list = [a3m_lines.split("\n") for a3m_lines in a3m_lines_mmseqs2]
    seq_list = []
    # Skip the first sequence, is as already been written
    for j in range(2, len(a3m_lines_list[0])):
        if a3m_lines_list[0][j].startswith(">"):
            title = ">"
            for i in range(len(a3m_lines_list)):
                title += " " + a3m_lines_list[i][j][1:]
        elif a3m_lines_list[0][j]:
            seq = ""
            for i in range(len(a3m_lines_list)):
                seq += "".join([aa for aa in a3m_lines_list[i][j] if aa in AA_CLASSIC])[
                    start_indice[i] : end_indice[i]
                ]
            if seq not in seq_list:
                f_out.write(f"{title}\n{seq}\n")
                seq_list.append(seq)

    print(f"- keep {len(seq_list)} sequences in the paired alignement")


def write_unpaired_alignement(
    f_out, a3m_lines_mmseqs2, start_indice, end_indice, chains_len
):
    """Write the unpaired alignement.

    Parameters
    ----------
    f_out : file
        The output file.
    a3m_lines_mmseqs2 : list of str
        The list of sequences.
    start_indice : list of int
        The list of start indices.
    end_indice : list of int
        The list of end indices.

    Returns
    -------
    None
    """

    # Unpaired alignement
    for i, a3m_lines in enumerate(a3m_lines_mmseqs2):
        line_list = a3m_lines.split("\n")

        insert_before = sum(chains_len[:i]) * "-"
        insert_after = sum(chains_len[i + 1 :]) * "-"

        seq_list = []

        for line in line_list:
            if line.startswith(">"):
                title = line
            elif line:
                seq = "".join([aa for aa in line if aa in AA_CLASSIC])[
                    start_indice[i] : end_indice[i]
                ]
                if not only_insert(seq) and seq not in seq_list:
                    f_out.write(f"{title}\n{insert_before}{seq}{insert_after}\n")
                    seq_list.append(seq)

                assert len(f"{insert_before}{seq}{insert_after}") == sum(
                    chains_len
                ), f"{seq} {lower_num} {sum(chains_len)} {len(f'{insert_before}{seq}{insert_after}')}"

        print(
            f"seq {i:2} : keep {len(seq_list):5} sequences in the unpaired alignement"
        )

def insert_seq(seq_ori, seq_new, gap_cost=-15):
    assert len(seq_new) > len(seq_ori)
    align = alignement.align_seq_cython(seq_ori, seq_new, gap_cost) 

    for insert_num in range(len(align[0])):
        if align[0][insert_num] != "-":
            break
    print(insert_num)
    
    seq_new = seq_new[:insert_num]+seq_ori+seq_new[len(seq_ori)+insert_num:]
    return seq_new

def create_full_alignement(
    name,
    original_seq_list,
    new_seq_list,
    out_dir,
    host_url="https://api.colabfold.com",
):
    """Create a full msa alignement from a list of sequences.

    Parameters
    ----------
    name : str
        The name of the alignement.
    original_seq_list : list of str
        The list of original sequences.
    new_seq_list : list of str
        The list of new sequences.
    out_dir : str
        The output directory.
    host_url : str
        The url of the host.

    Returns
    -------
    None
    """

    version = importlib_metadata.version("colabfold")
    user_agent = f"colabfold/{version}"

    # remove duplicates before searching
    ori_query_seqs_unique = []
    for seq in original_seq_list:
        if seq not in ori_query_seqs_unique:
            ori_query_seqs_unique.append(seq)
    new_query_seqs_unique = []
    for seq in new_seq_list:
        if seq not in new_query_seqs_unique:
            new_query_seqs_unique.append(seq)
    assert len(ori_query_seqs_unique) == len(
        new_query_seqs_unique
    ), f"len(ori_query_seqs_unique) = {len(ori_query_seqs_unique)} != len(new_query_seqs_unique) = {len(new_query_seqs_unique)}"

    chains_len = [len(seq) for seq in ori_query_seqs_unique]

    # determine how many times is each sequence is used (cardinality)
    query_seqs_cardinality = [0] * len(ori_query_seqs_unique)
    for seq in original_seq_list:
        seq_idx = ori_query_seqs_unique.index(seq)
        query_seqs_cardinality[seq_idx] += 1

    new_query_seqs_cardinality = [0] * len(ori_query_seqs_unique)
    for seq in new_seq_list:
        seq_idx = new_query_seqs_unique.index(seq)
        new_query_seqs_cardinality[seq_idx] += 1

    assert (
        query_seqs_cardinality == new_query_seqs_cardinality
    ), f"query_seqs_cardinality = {query_seqs_cardinality} != new_query_seqs_cardinality = {new_query_seqs_cardinality}"
    print(f"query_seqs_cardinality = {query_seqs_cardinality}")

    # check that original_seq is in new_seqs
    for ori_seq, new_seq in zip(ori_query_seqs_unique, new_query_seqs_unique):
        assert new_seq.find(ori_seq) != -1, f"{ori_seq} not in {new_seq}"

    result_dir = Path(out_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    print(f"Running MMseqs2 for {name} unpaired alignment")

    a3m_lines_mmseqs2_full = run_mmseqs2(
        new_query_seqs_unique,
        str(result_dir.joinpath(name + "_unpaired")),
        use_env=True,
        use_templates=False,
        use_pairing=False,
        host_url=host_url,
        user_agent=user_agent,
    )

    print(f"Running MMseqs2 for {name} paired alignment")

    a3m_lines_mmseqs2_paired = run_mmseqs2(
        new_query_seqs_unique,
        str(result_dir.joinpath(name + "_paired")),
        use_env=True,
        use_templates=False,
        use_pairing=True,
        host_url=host_url,
        user_agent=user_agent,
    )

    print("- Unpaired alignements")
    for i, a3m_line in enumerate(a3m_lines_mmseqs2_full):
        line_list = a3m_line.split("\n")
        print(
            f"seq {i} seq num = {len([line for line in line_list if line.startswith('>') ]):5},  res num ={len(line_list[1]):5}"
        )

    print("- Paired alignements")
    for a3m_line in a3m_lines_mmseqs2_paired:
        line_list = a3m_line.split("\n")
        print(
            f"seq {i} seq num = {len([line for line in line_list if line.startswith('>') ]):5},  res num ={len(line_list[1]):5}"
        )

    start_indice = []
    end_indice = []

    for i, a3m_lines in enumerate(a3m_lines_mmseqs2_full):
        line_list = a3m_lines.split("\n")
        start_indice.append(line_list[1].find(ori_query_seqs_unique[i]))
        end_indice.append(start_indice[-1] + len(ori_query_seqs_unique[i]))

    out_name = str(result_dir.joinpath(name + ".a3m"))
    with open(out_name, "w") as f_out:

        f_out.write(
            f'#{",".join([str(i) for i in chains_len])}	{",".join([str(i) for i in query_seqs_cardinality])}\n'
        )
        f_out.write(f">101	102\n")
        f_out.write(f'{"".join(ori_query_seqs_unique)}\n')

        write_paired_alignement(
            f_out, a3m_lines_mmseqs2_paired, start_indice, end_indice
        )
        write_unpaired_alignement(
            f_out, a3m_lines_mmseqs2_full, start_indice, end_indice, chains_len
        )


if __name__ == "__main__":

    original_seq = "GPEKEWVEQDEPGVYITLTALAGGARDLKRVRFSRKRFSEIQAEQWWADNRGRVYEQYNVRMV"
    # original_seq += ":GPEKEWVEQDEPGVYITLTALAGGARDLKRVRFSRKRFSEIQAEQWWADNRGRVYEQYNVRMV"
    original_seq += ":GPKWVKTDSDFIVLEI"
    original_seq_list = original_seq.split(":")
    # chains_len = [len(seq) for seq in original_seq_list]

    # new_seq = "GPEKEWVEQDEPGVYITLTALAGGARDLKRVRFSRKRFSEIQAEQWWADNRGRVYEQYNVRMV:MKFFNWMQNKLGGKQENRKSNTSTSTTYAKPEPREEFSDWPHSLLAIGTFGNNNEITQNIENQNTQQEDPSSSEEVPDFTPEEIGKLQKELTRLLRRKPNVEKEISELPLDRFLNCPSSLEVDRRISNALCSESGGDKDEDIEKTLSVILDKCKDICAEKSKKSIGKKSISFLLKKMFVCRSGFAPTPSLRDTLQESRMEKLLRTMLHKKLYTQNNSRAPVLKKCLENKKSIKKRNEDEAEERIDEGPKWVKTDSDFIVLEI"
    new_seq = (
        "MADLVTYSNADHNLEQALITLKKGTQLLKYGRKGKPKFYPFRLSSDEKSLIWISSSGEKRLKLASVSKIV"
        "PGQRTAVFQRYLRPEKDYLSFSLLYNGKKKSLDLICKDKVEAEIWIGGLKTLISTGQGGRSKIDGWSGGG"
        "LSVDASRELTSSSPSSSSASASRGHSSPGTPFNIDPITSPKSAEPEVPPTDSEKSHVALDNKNMQTKVSG"
        "SDGFRVSVSSAQSSSSHGSAADDSDALGDVYIWGEVICDNVVKVGIDKNASYLTTRTDVLVPKPLESNIV"
        "LDVHQIACGVRHAAFVTRQGEIFTWGEESGGRLGHGIGKDVFHPRLVESLTATSSVDFVACGEFHTCAVT"
        "LAGELYTWGDGTHNVGLLGHGSDISHWIPKRIAGSLEGLHVASVSCGPWHTALITSYGRLFTFGDGTFGV"
        "LGHGDKETVQYPREVESLSGLRTIAVSCGVWHTAAVVEIIVTQSNSSSVSSGKLFTWGDGDKNRLGHGDK"
        "DPRLKPTCVPALIDYNFHKIACGHSLTVGLTTSGQVFTMGSTVYGQLGNLQTDGKLPCLVEDKLASEFVE"
        "EISCGAYHVAALTSRNEVYTWGKGANGRLGHGDLEDRKVPTIVEALKDRHVKYIACGSNYTAAICLHKWV"
        "SGAEQSQCSTCRLAFGFTRKRHNCYNCGLVHCHSCSSKKAFRAALAPSAGRLYRVCDSCYVKLSKVSEIN"
        "DTNRRNSAVPRLSGENRDRLDKSEIRLAKFGTSNMDLIKQLDSKAAKQGKKTDTFSLGRNSQLPSLLQLK"
        "DAVQSNIGDMRRATPKLAQAPSGISSRSVSPFSRRSSPPRSATPMPSTSGLYFPVGIADNMKKTNEILNQ"
        "EIVKLRTQVDSLTQKCEFQEVELQNSVKKTQEALALAEEESAKSRAAKEAIKSLIAQLKDVAEKLPPGES"
        "VKLACLQNGLDQNGFHFPEENGFHPSRSESMTSSISSVAPFDFAFANASWSNLQSPKQTPRASERNSNAY"
        "PADPRLSSSGSVISERIEPFQFQNNSDNGSSQTGVNNTNGPEKEWVEQDEPGVYITLTALAGGARDLKRV"
        "RFSRKRFSEIQAEQWWADNRGRVYEQYNVRMVEKSTASQTHRDRDEEEEDIPH"
    )
    # new_seq += ":MADLVTYSNADHNLEQALITLKKGTQLLKYGRKGKPKFYPFRLSSDEKSLIWISSSGEKRLKLASVSKIV"\
    #          "PGQRTAVFQRYLRPEKDYLSFSLLYNGKKKSLDLICKDKVEAEIWIGGLKTLISTGQGGRSKIDGWSGGG"\
    #          "LSVDASRELTSSSPSSSSASASRGHSSPGTPFNIDPITSPKSAEPEVPPTDSEKSHVALDNKNMQTKVSG"\
    #          "SDGFRVSVSSAQSSSSHGSAADDSDALGDVYIWGEVICDNVVKVGIDKNASYLTTRTDVLVPKPLESNIV"\
    #          "LDVHQIACGVRHAAFVTRQGEIFTWGEESGGRLGHGIGKDVFHPRLVESLTATSSVDFVACGEFHTCAVT"\
    #          "LAGELYTWGDGTHNVGLLGHGSDISHWIPKRIAGSLEGLHVASVSCGPWHTALITSYGRLFTFGDGTFGV"\
    #          "LGHGDKETVQYPREVESLSGLRTIAVSCGVWHTAAVVEIIVTQSNSSSVSSGKLFTWGDGDKNRLGHGDK"\
    #          "DPRLKPTCVPALIDYNFHKIACGHSLTVGLTTSGQVFTMGSTVYGQLGNLQTDGKLPCLVEDKLASEFVE"\
    #          "EISCGAYHVAALTSRNEVYTWGKGANGRLGHGDLEDRKVPTIVEALKDRHVKYIACGSNYTAAICLHKWV"\
    #          "SGAEQSQCSTCRLAFGFTRKRHNCYNCGLVHCHSCSSKKAFRAALAPSAGRLYRVCDSCYVKLSKVSEIN"\
    #          "DTNRRNSAVPRLSGENRDRLDKSEIRLAKFGTSNMDLIKQLDSKAAKQGKKTDTFSLGRNSQLPSLLQLK"\
    #          "DAVQSNIGDMRRATPKLAQAPSGISSRSVSPFSRRSSPPRSATPMPSTSGLYFPVGIADNMKKTNEILNQ"\
    #          "EIVKLRTQVDSLTQKCEFQEVELQNSVKKTQEALALAEEESAKSRAAKEAIKSLIAQLKDVAEKLPPGES"\
    #          "VKLACLQNGLDQNGFHFPEENGFHPSRSESMTSSISSVAPFDFAFANASWSNLQSPKQTPRASERNSNAY"\
    #          "PADPRLSSSGSVISERIEPFQFQNNSDNGSSQTGVNNTNGPEKEWVEQDEPGVYITLTALAGGARDLKRV"\
    #          "RFSRKRFSEIQAEQWWADNRGRVYEQYNVRMVEKSTASQTHRDRDEEEEDIPH"
    new_seq += (
        ":MKFFNWMQNKLGGKQENRKSNTSTSTTYAKPEPREEFSDWPHSLLAIGTFGNNNEITQNIENQNTQQE"
        "DPSSSEEVPDFTPEEIGKLQKELTRLLRRKPNVEKEISELPLDRFLNCPSSLEVDRRISNALCSESGGDK"
        "DEDIEKTLSVILDKCKDICAEKSKKSIGKKSISFLLKKMFVCRSGFAPTPSLRDTLQESRMEKLLRTMLH"
        "KKLYTQNNSRAPVLKKCLENKKSIKKRNEDEAEERIDEGPKWVKTDSDFIVLEI"
    )

    new_seq_list = new_seq.split(":")

    create_full_alignement(
        "6LV0",
        original_seq_list,
        new_seq_list,
        out_dir="tmp",
    )

    # blastcmd IDENTIFIANT
    # chaque query -> evalue

    # blast -> pour un ensemble de donnée (nr)
