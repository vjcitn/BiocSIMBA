#' configuration list for training
#' @export
pbg_train_config = function(entity_path = "", edge_paths = "", checkpoint_path = "",
    entities = structure(list(), names = character(0)), relations = list(),
    dynamic_relations = FALSE, dimension = 50L, global_emb = FALSE,
    comparator = "dot", num_epochs = 10L, workers = 4L, num_batch_negs = 50L,
    num_uniform_negs = 50L, loss_fn = "softmax", lr = 0.1,
    early_stopping = FALSE, regularization_coef = 0, wd = 0,
    wd_interval = 50L, eval_fraction = 0.05, eval_num_batch_negs = 50L,
    eval_num_uniform_negs = 50L, checkpoint_preservation_interval = NULL) {

list(entity_path = entity_path, edge_paths = edge_paths, checkpoint_path = checkpoint_path, 
    entities = entities, relations = relations,
    dynamic_relations = dynamic_relations, dimension = dimension, global_emb = global_emb, 
    comparator = comparator, num_epochs = num_epochs, workers = workers, num_batch_negs = num_batch_negs, 
    num_uniform_negs = num_uniform_negs, loss_fn = loss_fn, lr = lr, 
    early_stopping = early_stopping, regularization_coef = regularization_coef, wd =wd, 
    wd_interval = wd_interval, eval_fraction = eval_fraction, eval_num_batch_negs = eval_num_batch_negs, 
    eval_num_uniform_negs = eval_num_uniform_negs, checkpoint_preservation_interval = checkpoint_preservation_interval)
}



#' after graph produced, train for embedding using PBG
#' @param output character(1)
#' @param auto_wd logical(1)
#' @param save_wd logical(1)
#' @param simconf like result to call to `simba_config()`
#' @param pbg_config list(), eventually to be handled as dict
#' @examples
#' h5path = "/home/vincent/tenx3k.h5ad"
#' g1 = gen_graph_sce(h5path)
#' tt = pbg_train_embedding( simconf=simba_config(workdir=g1$workdir) )
#' @export
pbg_train_embedding = function( output='model', auto_wd = TRUE,
    save_wd = TRUE,  simconf=simba_config(), pbg_config = pbg_train_config() ) {
  proc = basilisk::basiliskStart(bsklenv, testload="simba") # avoid package-specific import
  basilisk::basiliskRun(proc, function( output, auto_wd, save_wd, simconf, pbg_config) {
     sconf = simconf
     pconf = pbg_config
     sref = reticulate::import("simba")
     sref$settings$set_workdir(sconf$workdir)
#     sref$settings$pbg_params = dict(pconf)  # should pass
     sref$tl$pbg_train( auto_wd = auto_wd, save_wd = save_wd,
           output=output)
     list(simba_conf = sconf, pbg_conf = pconf )
     }, output, auto_wd, save_wd, simconf, pbg_config )
}

