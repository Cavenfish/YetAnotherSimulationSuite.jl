
struct Constraint{I,A} <: MyConstraint
  inds::I
  apply!::A
end

fixAtoms(atoms) = Constraint(atoms, fixAtoms!)

function fixAtoms!(forces, atoms)
  for i in atoms
    forces[i] .= 0.0
  end
end


