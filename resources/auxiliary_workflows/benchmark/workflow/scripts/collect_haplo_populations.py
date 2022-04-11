import os
import shutil

def main(list_dname_work, dname_out):

    if os.path.exists(dname_out)==False:
            os.makedirs(dname_out)

    dname_out = dname_out + "/"

    file1 = "/local_haplo_pairwise_distance.pdf"
    file2 = "/pairwise_hamming_distance.csv"

    for dname_work in list_dname_work:
        out_dir = dname_work.split('simulated_reads/')[1].split('/work')[0]
        out_path = dname_out + out_dir

        os.makedirs(out_path)

        src = dname_work + file1
        dst = out_path + file1
        shutil.copy(src, dst)  # dst can be a folder; use shutil.copy2() to preserve timestamp

        src = dname_work + file2
        dst = out_path + file2
        shutil.copy(src, dst)


if __name__ == "__main__":
    main(
        snakemake.input.dname_work,
        snakemake.output.dname_out,
    )
