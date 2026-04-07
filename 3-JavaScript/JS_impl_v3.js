import fs from "node:fs";
import { SimulationConfig } from "./JS_v3_config.js";
import {
  BoundaryConditions,
  CH_Solver,
  DoubleWellBulk,
  GradientEnergy,
  Mesh,
  RandomInit,
} from "./JS_v3_data_structures.js";

const configPath = process.argv[2] ?? "simu_config.json";
const simulationConfig = SimulationConfig.fromFile(configPath);
const boundaryConditions = new BoundaryConditions(
  simulationConfig.mesh.boundaryConditions,
);

fs.mkdirSync(simulationConfig.output.folder_name, { recursive: true });
console.log(`Using config: ${configPath}`);
console.log(`Output folder: ${simulationConfig.output.folder_name}`);

const con_mesh = new Mesh(
  simulationConfig.mesh.Nx,
  simulationConfig.mesh.Ny,
  simulationConfig.mesh.dx,
  simulationConfig.mesh.dy,
  undefined,
  boundaryConditions,
);
con_mesh.init_geometry(
  new RandomInit(
    simulationConfig.mesh.con_init,
    -simulationConfig.mesh.dcon,
    simulationConfig.mesh.dcon,
  ),
);
con_mesh.boundCons.applyAll(con_mesh);

const bulk = new DoubleWellBulk(simulationConfig.energy.A);
const interfacial = new GradientEnergy(simulationConfig.energy.kappa);

const solver = new CH_Solver(
  con_mesh,
  simulationConfig.solver.M,
  bulk.dbulk_dc(),
  interfacial.dint_dc(),
  undefined,
  simulationConfig.solver.dt,
);

function assertNoNaN(mesh, step) {
  const firstNaN = mesh.findFirstNaN();
  if (firstNaN !== undefined) {
    throw new Error(
      `NaN detected at step ${step} at position (x=${firstNaN.x}, y=${firstNaN.y}), flat index ${firstNaN.idx}.`,
    );
  }
}

const nanCheckEvery = simulationConfig.output.nanCheckEvery;
if (nanCheckEvery !== false) {
  assertNoNaN(con_mesh, -1);
}

performance.mark("start-simulation");
for (let istep = 0; istep < simulationConfig.solver.Nstep + 1; istep++) {
  solver.step();
  if (nanCheckEvery !== false && istep % nanCheckEvery === 0) {
    assertNoNaN(con_mesh, istep);
    console.log(`NaN check passed at step ${istep}.`);
  }

  if (istep % simulationConfig.output.writeEvery === 0) {
    con_mesh.write_VTK(istep, simulationConfig.output.folder_name);
  }
  if (istep % simulationConfig.output.logEvery === 0) {
    console.log(`Step: ${istep}`);
  }
}
performance.mark("end-simulation");
performance.measure(
  "Total Simulation Time",
  "start-simulation",
  "end-simulation",
);

console.log(
  `Total Simulation Time: ${performance.getEntriesByName("Total Simulation Time")[0].duration.toFixed(2)} ms`,
);
