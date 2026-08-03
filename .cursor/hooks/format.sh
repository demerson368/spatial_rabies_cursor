input=$(cat)
file_path=$(echo "$input" | jq -r '.file_path // empty')
if [[ "$file_path" == *.ts || "$file_path" == *.tsx ]]; then
  npx prettier --write "$file_path" 2>/dev/null
fi
echo '{}'
exit 0

