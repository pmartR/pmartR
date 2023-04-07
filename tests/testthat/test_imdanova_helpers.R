context("test imd anova helper fns")

test_that("comparison matrix is properly constructed", {
    base_name = sample(LETTERS, 1)
    gname1 = paste0(base_name, sample(LETTERS, 5), collapse = "")
    gname2 = paste(gname1, paste(sample(LETTERS, 5), collapse=""), sep = "")
    gvec1 = rep(gname1, 5)
    gvec2 = rep(gname2, 5)

    group_vec = c(
        rep(base_name, 5),
        gvec1,
        gvec2,
        rep(base_name, 5),
        gvec1,
        gvec2
    )

    group_df = structure(
        list(SampleID = make.unique(sample(LETTERS, 30, replace = TRUE)), 
            Group = group_vec), 
        main_effects = "Group", 
        nonsingleton_groups = c(gvec1[1], gvec2[1], base_name), 
        row.names = c(NA, 30L), class = "data.frame")

    comps1 = structure(list(Control = c(base_name, base_name, gname1), Test = c(gname1, 
    gname2, gname2)), class = "data.frame", row.names = c(NA, -3L))

    comps2 = structure(list(Control = c(gname1, base_name, gname1), Test = c(gname2,
    gname2, base_name)), class = "data.frame", row.names = c(NA, -3L))

    comps3 = structure(list(Control = c(gname2, gname2, base_name), Test = c(gname1,
    base_name, gname1)), class = "data.frame", row.names = c(NA, -3L))

    mat1 = create_c_matrix(group_df, comps1)
    mat2 = create_c_matrix(group_df, comps2)
    mat3 = create_c_matrix(group_df, comps3)

    res1 = list(
        cmat = matrix(c(-1, 1, 0, -1, 0, 1, 0, -1, 1),nrow=3, byrow=TRUE),
        names = c(
            paste0(gname1, "_vs_", base_name),
            paste0(gname2, "_vs_", base_name),
            paste0(gname2, "_vs_", gname1)
        )
    )

    res2 = list(
        cmat = matrix(c(0, -1, 1, -1, 0, 1, 1, -1, 0),nrow=3, byrow=TRUE),
        names = c(
            paste0(gname2, "_vs_", gname1),
            paste0(gname2, "_vs_", base_name),
            paste0(base_name, "_vs_", gname1)
        )
    )

    res3 = list(
        cmat = matrix(c(0, 1, -1, 1, 0, -1, -1, 1, 0),nrow=3, byrow=TRUE),
        names = c(
            paste0(gname1, "_vs_", gname2),
            paste0(base_name, "_vs_", gname2),
            paste0(gname1, "_vs_", base_name)
        )
    )

    expect_equal(mat1, res1)
    expect_equal(mat2, res2)
    expect_equal(mat3, res3)
})

