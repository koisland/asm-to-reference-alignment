def read_sample_tbl():
    df = pd.read_csv(config.get("tbl"), sep="\t")
    df.asm = df.asm.map(os.path.abspath)
    df["asm"] = df.asm.str.split(",")
    df = df.explode("asm")
    df["num"] = df.groupby(level=0).cumcount() + 1
    if len(df["num"].unique()) > 1:
        new_index = df["sample"] + "_" + df["num"].astype(str)
    else:
        new_index = df["sample"].astype(str)

    df.set_index(new_index, inplace=True)
    return df


def get_asm(wc):
    return df.loc[str(wc.sm)].asm


def get_ref(wc):
    return config.get("ref")[wc.ref]


def get_fai(wc):
    return config.get("ref")[wc.ref] + ".fai"
